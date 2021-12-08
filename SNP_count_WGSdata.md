# Counting sex-limited and sex-shared SNPs in WGS data from west

```
#!/usr/bin/env perl
use strict;
use warnings;

# 	This program executes another perl script called
#   XB_SNPcount_EastWest_limits.pl with consecutive
#   windows and pastes output into the same file

my $status;
my $end;

my $x;
for ($x = 0 ; $x <= 57000000 ; $x+=100000) {
	$end = $x+99999;
	$status = system("./XB_SNPcount_EastWest_limits.pl allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz.recode.vcf.tab 5_6_1_2_4_3 EastWest_windows_counts.txt ".$x." ".$end);
}
```

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

# to execute type XB_SNPcount_EastWest_limits.pl inputfile.tab 5_6_1_2_4_3 output_counts.txt lower upper

# where the 5_6_1_2_4_3 indicates that the 3rd and 4th individuals are the East_female (1)
# and East_male (2), the sixth and fifth individuals are the Lab_female (3) and the Lab_male (4),
# and the first and second individual are the West_female (5), and West_male (6).

# This is useful because it enables me to switch around the counts for comparison/control (e.g. I can
# count how many male-specific SNPs from the East and Lab are in the west)

# lower and upper are the coordinates that we want to use.  For the whole SL region, use 0 and 54000000

my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $outputfile = $ARGV[2];
my $lower = $ARGV[3];
my $upper = $ARGV[4];

#my $SL_limit = 54000000; # This is the upper limit of the SL region; 54million for Chr8L for XB
					  # Lower for Chr7S 

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file 1.\n";
	exit;
}


my @individuals=split("_",$input2); # This will tell us the order of individuals
my $a;


unless (open(OUTFILE, ">>$outputfile"))  {
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
my $West_female_heteroz_and_West_male_sex_shared_specific_homoz=0;
my $West_male_heteroz_and_West_female_sex_shared_specific_homoz=0;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if(($temp[0] ne '#CHROM')&&($temp[1] > $lower)&&($temp[1] < $upper)){
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
						# Now count the number of sites where the west female is heteroz and the west male is homoz for the sex-shared SNP
						if(($West_female_allele_1 ne $West_female_allele_2)&&(($West_male_allele_1 eq $sex_shared_SNP)&&($West_male_allele_2 eq $sex_shared_SNP))){
							$West_female_heteroz_and_West_male_sex_shared_specific_homoz+=1;
						}
						elsif(($West_male_allele_1 ne $West_male_allele_2)&&(($West_female_allele_1 eq $sex_shared_SNP)&&($West_female_allele_2 eq $sex_shared_SNP))){
							$West_male_heteroz_and_West_female_sex_shared_specific_homoz+=1;
						}
					} # end of check for het and homoz East and Lab
				} # end of check for nucleotides
			} # end of check that position is biallelic with one character alleles
		 } # end of check for 14 columns in @temp after split	
	} # end of check for beginning of the line or beyond SL limit
} # end of while statement	

#print @individuals,"\n";

# Print results to an output file
print $lower," ",$upper,"\n"; #"Lower\tUpper\tIndividual\tn_homoz_female_specific\tn_homoz_sex_shared\tn_het\tIndividual\tn_homoz_female_specific\tn_homoz_sex_shared\tn_het\tConserved\tWest_female_heteroz_and_West_male_sex_shared_specific_homoz\tControl\tWest_male_heteroz_and_West_female_sex_shared_specific_homoz\n";
print OUTFILE $lower,"\t",$upper,"\tWest_female\t", $West_female_female_specific_homoz,"\t", $West_female_sex_shared_specific_homoz,"\t", $West_female_heteroz,"\t";
print OUTFILE "West_male\t", $West_male_female_specific_homoz,"\t", $West_male_sex_shared_specific_homoz,"\t", $West_male_heteroz,"\tConserved\t",$West_female_heteroz_and_West_male_sex_shared_specific_homoz,"\tControl\t",$West_male_heteroz_and_West_female_sex_shared_specific_homoz,"\n";
close OUTFILE;
```
