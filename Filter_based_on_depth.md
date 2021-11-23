# Filtering based on rolling average of depth

Get depth
```
#!/bin/sh
#SBATCH --job-name=VariantsToTable
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=10gb
#SBATCH --output=VariantToTable.%J.out
#SBATCH --error=VariantToTable.%J.err
#SBATCH --account=def-ben


# This script will execute the GATK command "VariantsToTable

# execute like this:
# sbatch 2021_VariantFiltration.sh path_and_file

module load nixpkgs/16.09 gatk/4.1.0.0

gatk --java-options -Xmx8G VariantsToTable -V ${1} -F CHROM -F POS -GF DP -O ${1}.table
```
Load R
```
module load StdEnv/2020 java/13.0.2 r/4.1.2
```
Need to change file names twice in this script
```
library (ggplot2)
setwd("./")
# get all the files with the site data
files <- list.files(path = ".", pattern = "allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz.table")

my_data<-c()
temp<-c()

# read in the data and name the df based on the file name
for(f in 1:length(files)) {
  temp <- read.table(files[f], header = T)
  my_data <- rbind(my_data,temp)
  temp<-c()
}  

# my_data now has coverage per site info for each sample for all sites, genomewide
dim(my_data)
colnamez <- colnames(my_data)

# here is a function to calculate moving averages
# https://stackoverflow.com/questions/743812/calculating-moving-average
moving_fun <- function(x, w, FUN, ...) {
  # x: a double vector
  # w: the length of the window, i.e., the section of the vector selected to apply FUN
  # FUN: a function that takes a vector and return a summarize value, e.g., mean, sum, etc.
  # Given a double type vector apply a FUN over a moving window from left to the right, 
  #    when a window boundary is not a legal section, i.e. lower_bound and i (upper bound) 
  #    are not contained in the length of the vector, return a NA_real_
  if (w < 1) {
    stop("The length of the window 'w' must be greater than 0")
  }
  output <- x
  for (i in 1:length(x)) {
    # plus 1 because the index is inclusive with the upper_bound 'i'
    lower_bound <- i - w + 1
    if (lower_bound < 1) {
      output[i] <- NA_real_
    } else {
      output[i] <- FUN(x[lower_bound:i, ...])
    }
  }
  output
}

# example
# v <- seq(1:10)

# compute a MA(2)
# moving_fun(v, 2, mean)


# make a new dataframe that has the moving average for each sample throughout the 
# genome.  No worries if the window goes across chromosomes
mv_ave_df2 <- c()
for (i in 3:length(my_data)){  
  # calculate the moving average for each column
  moving_average <- moving_fun(my_data[,i], 50, mean)
  # add this to a new dataframe
  mv_ave_df2 <- cbind(mv_ave_df2,moving_average)
  # rename the column to match the sample
  colnames(mv_ave_df2)[i-2] = colnamez[i]
}  

# add chromosome and position data to mv_ave_df2
mv_ave_df3 <- data.frame(mv_ave_df2,my_data$CHROM,my_data$POS)

colnames(mv_ave_df3)[7] <- "CHROM"
colnames(mv_ave_df3)[8] <- "POS"


# calculate mean depth and sd per sample
mean_depth<- c()
sd_depth<- c()

for (i in 3:length(my_data)){
  x<- mean(my_data[,i],na.rm=T)
  mean_depth<-append(mean_depth,x, after = length(mean_depth))
  y<-sd(my_data[,i],na.rm=T)
  sd_depth<-append(sd_depth,y, after = length(sd_depth))
}  

cutoff_vector <- c()


# identify chr and positions in any sample that are >3 sd above that samples mean coverage
# based on the rolling average
# first make a vector with cutoff values for each sample
for (i in 3:ncol(mv_ave_df3)){
  cutoff <- mean_depth[i-2] + 3*sd_depth[i-2]
  cutoff_vector <- append(cutoff_vector,cutoff, after=length(cutoff_vector))
}

mean_over_cuttoff <- mean_depth/cutoff_vector
par(mfrow=c(2,2), mar=c(2,2,2,2)) 
plot(mean_over_cuttoff)
# give the elements in cutoff_vector some names
names(cutoff_vector) <- c('BJE4441_female_west_merge_sorted_dedup.bam.DP',
                          'BJE4442_male_west_merge_sorted_dedup.bam.DP',
                          'BJE4515_female_east_merge_sorted_dedup.bam.DP',
                          'BJE4536_male_east_merge_sorted_dedup.bam.DP',
                          'SRR6357672.1_trim_fixed_merged_dedup.bam.DP',
                          'SRR6357673_trim_sorted_dedup.bam.DP' )

# now cycle through each column and identify chr and pos of the bad ones
sub <- NULL
sub <- subset(mv_ave_df3, 
              BJE4441_female_west_merge_sorted_dedup.bam.DP > cutoff_vector["BJE4441_female_west_merge_sorted_dedup.bam.DP"]  | 
                BJE4442_male_west_merge_sorted_dedup.bam.DP > cutoff_vector["BJE4442_male_west_merge_sorted_dedup.bam.DP"] |
                BJE4515_female_east_merge_sorted_dedup.bam.DP  > cutoff_vector["BJE4515_female_east_merge_sorted_dedup.bam.DP"] |
                BJE4536_male_east_merge_sorted_dedup.bam.DP  > cutoff_vector["BJE4536_male_east_merge_sorted_dedup.bam.DP"] |
                SRR6357672.1_trim_fixed_merged_dedup.bam.DP  > cutoff_vector["SRR6357672.1_trim_fixed_merged_dedup.bam.DP"] |
                SRR6357673_trim_sorted_dedup.bam.DP > cutoff_vector["hecki_PF647.DP"] |
                SRR6357673_trim_sorted_dedup.bam.DP> cutoff_vector["SRR6357673_trim_sorted_dedup.bam.DP"] 
              )

dim(sub)
dim(mv_ave_df3)


to_file <- data.frame(sub$CHROM, sub$POS)

# write to file
write.table(to_file, "Chr8L_positions_to_exclude.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = F)
```

Remove quotes
```
sed -i 's/"//g' Chr8L_positions_to_exclude.txt
```
Now remove these sites
```
module load nixpkgs/16.09 intel/2018.3 vcftools/0.1.16
vcftools --gzvcf allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz --out allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz.recode.vcf.gz --exclude-positions Chr8L_positions_to_exclude.txt --recode
```

make a tab file
```
module load nixpkgs/16.09  intel/2018.3 vcftools/0.1.16
vcf-to-tab < in.vcf > out.tab
```
Count SNPs
```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use List::Util qw( min max );


#  This program reads in a tab delimited genotype file generated
#  from vcftools that has a West_female, West_male, East_female,
# East_male, Lab_male_72, Lab_female_73
#  in that order.

# The script will search for binucleotide sites with no missing data

# It will then find positions where the East_female and the Lab_female
# are heterozygous but the East_male and the Lab_male are homozygous

# Then it will count how frequently the female-specific SNP occurs in 
# the West_female and the West_male and also how frequently the sex-shared
# SNP occurs in the West_female and the West_male 

# Data should be from one chromosome.

# on graham, you need to load perl first:
# module load StdEnv/2020 perl/5.30.2

# to execute type XB_SNPcountEastWest.pl inputfile.tab 5_6_1_2_4_3 output_counts.txt 

# where the 5_6_1_2_4_3 indicates that the 3rd and 4th individuals are the East_female (1)
# and East_male (2), the sixth and fifth individuals are the Lab_female (3) and the Lab_male (4),
# and the first and second individual are the West_female (5), and West_male (6).

# This is useful because it enables me to switch around the counts for comparison/control (e.g. I can
# count how many male-specific SNPs from the East and Lab are in the west)


my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $outputfile = $ARGV[2];

my $SL_limit = 54000000; # This is the upper limit of the SL region; 54million for Chr8L for XB
					  # Lower for Chr7S 

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file 1.\n";
	exit;
}


my @individuals=split("_",$input2); # This will tell us the order of individuals
my $a;


unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

my $allele1;
my $allele2;
my $allele3;
my $allele4;
my $allele5;
my $allele6;
my $allele7;
my $allele8;
my $allele9;
my $allele10;
my $allele11;
my $allele12;

my $counter=0;
my $concat="";

my $East_female_allele_1;
my $East_female_allele_2;
my $East_male_allele_1;
my $East_male_allele_2;

my $Lab_female_allele_1;
my $Lab_female_allele_2;
my $Lab_male_allele_1;
my $Lab_male_allele_2;

my $West_female_allele_1;
my $West_female_allele_2;
my $West_male_allele_1;
my $West_male_allele_2;

my @temp;
my @all_alleles;
my @all_alleles_uniq;
my $female_specific_SNP;
my $sex_shared_SNP;

# SNP counts
my $West_female_female_specific_homoz=0;
my $West_female_sex_shared_specific_homoz=0;
my $West_female_heteroz=0;
my $West_male_female_specific_homoz=0;
my $West_male_sex_shared_specific_homoz=0;
my $West_male_heteroz=0;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if(($temp[0] ne '#CHROM')&&($temp[1] < $SL_limit)){
		# add a check that there are only 14 columns after split
		if($#temp == 14){
		# load the alleles
			# load each of the 12 alleles
			# The identity of each allele needs to be assigned based on which
			# individual I have designated in the input string
			$counter=0; # reset the counter
			$concat=""; # reset the string for checking length
			for($a=3; $a<=$#temp; $a=$a+2){
				if($individuals[$counter] == 1){
					$East_female_allele_1 = $temp[$a];
					$East_female_allele_2 = $temp[$a+1];
					$concat=$concat.$temp[$a].$temp[$a+1];
				}
				elsif($individuals[$counter] == 2){
					$East_male_allele_1 = $temp[$a];
					$East_male_allele_2 = $temp[$a+1];
					$concat=$concat.$temp[$a].$temp[$a+1];
				}
				elsif($individuals[$counter] == 3){
					$Lab_female_allele_1 = $temp[$a];
					$Lab_female_allele_2 = $temp[$a+1];
					$concat=$concat.$temp[$a].$temp[$a+1];
				}
				elsif($individuals[$counter] == 4){
					$Lab_male_allele_1 = $temp[$a];
					$Lab_male_allele_2 = $temp[$a+1];
					$concat=$concat.$temp[$a].$temp[$a+1];
				}
				elsif($individuals[$counter] == 5){
					$West_female_allele_1 = $temp[$a];
					$West_female_allele_2 = $temp[$a+1];
					$concat=$concat.$temp[$a].$temp[$a+1];
				}
				elsif($individuals[$counter] == 6){
					$West_male_allele_1 = $temp[$a];
					$West_male_allele_2 = $temp[$a+1];
					$concat=$concat.$temp[$a].$temp[$a+1];
				}
				$counter+=1;
			}
			# check if there are only single positions and with only 2 nucleotide variants
			@all_alleles=split(//,$concat);
			@all_alleles_uniq = uniq(@all_alleles);
			if(($#all_alleles == 11)&&($#all_alleles_uniq == 1)){ # position is biallelic with one character alleles
				# now check that there are no gaps, asterisks, or other weird characters
				if(
					(($all_alleles_uniq[0] eq 'A')||($all_alleles_uniq[0] eq 'C')||($all_alleles_uniq[0] eq 'G')||($all_alleles_uniq[0] eq 'T')) &&
					(($all_alleles_uniq[1] eq 'A')||($all_alleles_uniq[1] eq 'C')||($all_alleles_uniq[1] eq 'G')||($all_alleles_uniq[1] eq 'T'))
					){ # this position is biallelic for two nucleotides
					# check if the East_female and Lab_female are heterozygous and the East_male and Lab_male are homozygous
					if(
					(($East_female_allele_1 ne $East_female_allele_2)&&($Lab_female_allele_1 ne $Lab_female_allele_2))&&
					(($East_male_allele_1 eq $East_male_allele_2)&&($Lab_male_allele_1 eq $Lab_male_allele_2)&&($East_male_allele_1 eq $Lab_male_allele_1))
					){	
						#print $line,"\n";
						# now identify the sex-shared and female-specific SNPs
						if($East_female_allele_1 eq $East_male_allele_1){
							$sex_shared_SNP = $East_female_allele_1;
							$female_specific_SNP = $East_female_allele_2;							
						} # end of if to assign female-specific and sex-shared SNP
						else{
							$sex_shared_SNP = $East_female_allele_2;
							$female_specific_SNP = $East_female_allele_1;
							# add a check for weirdness
							if($East_female_allele_2 ne $East_male_allele_1){
								print "Problem - the sex-shared SNP is not sex shared!\n";
							}
						} # end of else to assign female-specific and sex-shared SNP
					
						# Now count homoz and hets for each allele in the west female individuals
						if(($West_female_allele_1 eq $female_specific_SNP)&&($West_female_allele_2 eq $female_specific_SNP)){
								# The West_female is homoz for the female-specific SNP
								$West_female_female_specific_homoz+=1;
						}
						elsif(($West_female_allele_1 eq $sex_shared_SNP)&&($West_female_allele_2 eq $sex_shared_SNP)){
								# The West_female is homoz for the sex shared SNP
								$West_female_sex_shared_specific_homoz+=1;
								
						}
						elsif($West_female_allele_1 ne $West_female_allele_2){
								# The West_female is heterozygous for the female-specific and sex shared SNP
								$West_female_heteroz+=1;
						}
						else{
							print "Problem with assigning genotypes to West_female\n";
						}
						
						# Now count homoz and hets for each allele in the west male individuals
						if(($West_male_allele_1 eq $female_specific_SNP)&&($West_male_allele_2 eq $female_specific_SNP)){
								# The West_male is homoz for the female-specific SNP
								$West_male_female_specific_homoz+=1;								
						}
						elsif(($West_male_allele_1 eq $sex_shared_SNP)&&($West_male_allele_2 eq $sex_shared_SNP)){
								# The West_male is homoz for the sex shared SNP
								$West_male_sex_shared_specific_homoz+=1;
						}
						elsif($West_male_allele_1 ne $West_male_allele_2){
								# The West_male is heterozygous for the female-specific and sex shared SNP
								$West_male_heteroz+=1;
						}
						else{
							print "Problem with assigning genotypes to West_male\n";
						}
					} # end of check for het and homoz East and Lab
				} # end of check for nucleotides
			} # end of check that position is biallelic with one character alleles
		 } # end of check for 14 columns in @temp after split	
	} # end of check for beginning of the line or beyond SL limit
} # end of while statement	

print @individuals,"\n";

# Print results to stout
print "Individual\tn_homoz_female_specific\tn_homoz_sex_shared\tn_het\n";
print "West_female\t", $West_female_female_specific_homoz,"\t", $West_female_sex_shared_specific_homoz,"\t", $West_female_heteroz,"\n";
print "West_female\t", $West_male_female_specific_homoz,"\t", $West_male_sex_shared_specific_homoz,"\t", $West_male_heteroz,"\n";
# Print results to an output file
print OUTFILE "Individual\tn_homoz_female_specific\tn_homoz_sex_shared\tn_het\n";
print OUTFILE "West_female\t", $West_female_female_specific_homoz,"\t", $West_female_sex_shared_specific_homoz,"\t", $West_female_heteroz,"\n";
print OUTFILE "West_male\t", $West_male_female_specific_homoz,"\t", $West_male_sex_shared_specific_homoz,"\t", $West_male_heteroz,"\n";
close OUTFILE;
```
