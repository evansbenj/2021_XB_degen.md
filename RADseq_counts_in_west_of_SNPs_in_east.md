# RADseq counts in west of SNPs in east

Two expectations under the scenario of turnover in XB are that (1) the W chr should be fixed in the west and (2) there should be a new Y in the west.  I wrote a script to count female-specific SNPs in the east and then calculate their frequencies in females and males from the west.  It takes a tab delimited file as input.

First make a tab delimited file and subset it to include only the SL region:
```
module load StdEnv/2020 gcc/9.3.0 bcftools/1.11
bcftools view mpileup_raw_wildBorealis_AustinGenome.vcf.gz --regions Chr8L:1-54000000 -O z -o mpileup_raw_wildBorealis_AustinGenome_Chrs8_SL_only.vcf.gz

```

Now run this script:
```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use List::Util qw( min max );



#  This program reads in a tab delimited genotype file generated
#  from vcftools.  It will count SNPs that are specific to a focal population 
#  in this case east females, and the SNPs are putatively W-specific. These SNPs are not present
#  in the partners of the focal population, which in this case is east males
#  Then quantify how frequently the focal population SNPs are fixed in the each sex in the non-focal
#  population (west kenya M and Fs). Also calculate how frequently the putative Z SNPs are fixed in the west.

#  The prediction is that W-specific SNPs from the east will be fixed in the West because the Z is degenerate

#  Data should be from one chromosome, and probably wise to do it for the SL region and the non-SL region seperately

#  The prediction only applies for the SL region. For the non-SL region, a W-specific SNP from the east should
#  rarely be fixed in the west.

#  To deal with missing data and/or polymorphisms on the W, could quantify W-specific SNPs that are not fixed
#  in all east females or that have missing genotypes in some east females.  These should be more common than the
#  equivalent count in the non-Sex linked region of chr8L or elsewhere.

#  on cedar, you (may) need to load perl first:
#  module load StdEnv/2020 perl/5.30.2

#  to execute type Counts_east_SNPs_in_west.pl inputfile.tab output_counts_file positions_of_east_females positions_of_east_males positions_of_west_females positions_of_west_males 

# Counts_east_SNPs_in_west.pl temp.tab temp.out 27_28_29_30_31 52_53_54_55_56 3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18 32_33_34_35_36_37_38_39_40_41_42_43_44_45_46

# for Furman's XB file:
# ./Counts_east_SNPs_in_west.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L_sexlinked_only.tab temp.out 41_42_43_44_45 86_87_88_89_90 3_4_5_6_7_8_9_10_17_18_26_27_28_29 46_47_48_49_50_51_52_57_64_65_66_69_75_76_77

# for comparison between Nairobi/Nakuru to west:
# ./Counts_east_SNPs_in_west.pl mpileup_raw_justBorealis_allChrs.vcf_chr8L_sexlinked_only.tab temp.out 19_20_21_22 48 3_4_5_6_7_8_9_10_17_18_26_27_28_29 46_47_48_49_50_51_52_57_64_65_66_69_75_76_77


# where positions are the zero-indexed column number separated by underscores beginning with column 3 for the first sample
# Wundyani females: 27_28_29_30_31 (this is really columns 28_29_30_31_32 but we're using zero indexes)
# Wundyani males: 52_53_54_55_56

# extreme_west_females: 3_4_5_6_7_8_9_10_11_12_13_14_15_16_17_18
# extreme_west_males: 32_33_34_35_36_37_38_39_40_41_42_43_44_45_46

# central_females: 19_20_21_22_23_24_25_26
# central_males: 47_48_49_50_51

my $inputfile = $ARGV[0];
my $outputfile2 = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE2, ">$outputfile2"))  {
	print "I can\'t write to $outputfile2\n";
	exit;
}
print "Creating output file: $outputfile2\n";

my @east_female_indexes=split("_",$ARGV[2]);
my @east_male_indexes=split("_",$ARGV[3]);
my @west_female_indexes=split("_",$ARGV[4]);
my @west_male_indexes=split("_",$ARGV[5]);

my @east_female_nucleotides=();
my @east_male_nucleotides=();
my @west_female_nucleotides=();
my @west_male_nucleotides=();

my @uniq_east_female_nucleotides=();
my @uniq_east_male_nucleotides=();
my @uniq_west_female_nucleotides=();
my @uniq_west_male_nucleotides=();



# these arrays now have the zero-indexed column numbers of the four groups of genotypes

my @alleles;
my @temp;
my @lengths;
my $max;
my $a;
my $b;
my $start;
my $switch=0;

my $count_W_SNPs_in_east=0;
my $count_W_SNPs_in_east_fems=0;
my $count_W_SNPs_in_east_males=0;

my $count_Z_SNPs_in_east=0;
my $count_Z_SNPs_in_east_fems=0;
my $count_Z_SNPs_in_east_males=0;

my $count_W_SNPs_in_west=0; # this is the count of all W specific SNPs
my $count_W_SNPs_in_west_fems=0; # these are counters used for each position
my $count_W_SNPs_in_west_males=0; # these are counters used for each position

my $count_Z_SNPs_in_west=0;  # this is the count of all W specific SNPs in west
my $count_Z_SNPs_in_west_fems=0; # these are counters used for each position
my $count_Z_SNPs_in_west_males=0; # these are counters used for each position

my $W_specific;
my $Z_specific;
my $number_of_alleles;
my $number_of_alleles_f;
my $number_of_alleles_m;
my $frequency_of_east_W_SNPs_in_west;
my $frequency_of_east_Z_SNPs_in_west;
my $frequency_of_east_W_SNPs_in_west_fems;
my $frequency_of_east_W_SNPs_in_west_males;
my $frequency_of_east_Z_SNPs_in_west_fems;
my $frequency_of_east_Z_SNPs_in_west_males;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line); # split into individual nucleotides, even within a genotype
	if($temp[0] ne '#CHROM'){
			# load each of the alleles, the first three values are the chr, position, and ref
			@alleles=@temp;	
			# now see which is the longest
			@lengths=();
			for($a=3; $a<=$#alleles; $a++){
				push(@lengths,length($alleles[$a]));
			}
			$max = max @lengths;
			#print $max;
			# pad alleles if needed
			my @temp_alleles;
			# replace missing nucleotides with a dash
			for($a=3; $a<=$#alleles; $a++){
				if($alleles[$a] eq '.'){
					$alleles[$a] = "-";
				}	
			}	
			for($a=3; $a<=$#alleles; $a++){
				@temp_alleles=split(//,$alleles[$a]);
				if(($#temp_alleles+1) < $max){
					for($b=($#temp_alleles+1); $b<$max; $b++){
						$alleles[$a]=$alleles[$a]."-"; # add a gap to the end of short ones
					}	
				}	
			}		
			# now the lengths should all be the same
			
			# make arrays that have all the non-missing variation from each of the four populations
			@east_female_nucleotides=();
			@east_male_nucleotides=();
			@west_female_nucleotides=();
			@west_male_nucleotides=();
			
			for($a=0; $a<=$#east_female_indexes; $a++){
				if(($alleles[$east_female_indexes[$a]] ne ".")&&
				($alleles[$east_female_indexes[$a]] !~ /-/)&&
				($alleles[$east_female_indexes[$a]] ne "*")){
					push(@east_female_nucleotides,$alleles[$east_female_indexes[$a]]);
				}	
			}	
			for($a=0; $a<=$#east_male_indexes; $a++){
				if(($alleles[$east_male_indexes[$a]] ne ".")&&
				($alleles[$east_male_indexes[$a]] !~ /-/)&&
				($alleles[$east_male_indexes[$a]] ne "*")){
					push(@east_male_nucleotides,$alleles[$east_male_indexes[$a]]);
				}	
			}	
			for($a=0; $a<=$#west_female_indexes; $a++){
				if(($alleles[$west_female_indexes[$a]] ne ".")&&
				($alleles[$west_female_indexes[$a]] !~ /-/)&&
				($alleles[$west_female_indexes[$a]] ne "*")){
					push(@west_female_nucleotides,$alleles[$west_female_indexes[$a]]);
				}	
			}	
			for($a=0; $a<=$#west_male_indexes; $a++){
				if(($alleles[$west_male_indexes[$a]] ne ".")&&
				($alleles[$west_male_indexes[$a]] !~ /-/)&&
				($alleles[$west_male_indexes[$a]] ne "*")){
					push(@west_male_nucleotides,$alleles[$west_male_indexes[$a]]);
				}	
			}	
			
			# ok we now should have 4 arrays, each with non-missing nucleotides for each individual in each population
			# check the sizes of these arrays - usually should be even 
			print $#east_female_nucleotides+1," ",$#east_male_nucleotides+1," ",
				$#west_female_nucleotides+1," ",$#west_male_nucleotides+1," even\n";
			# check if there is polymorphism in the focal population (usually east females)
			@uniq_east_female_nucleotides = uniq @east_female_nucleotides;
			@uniq_east_male_nucleotides = uniq @east_male_nucleotides;
			@uniq_west_female_nucleotides = uniq @west_female_nucleotides;
			@uniq_west_male_nucleotides = uniq @west_male_nucleotides;
			
			# Now identify putative W-specific polymorphisms in the east females
			# first check if there is polymorphism in the east females after excluding missing data
			if(($#uniq_east_female_nucleotides == 1)&&($#uniq_east_male_nucleotides == 0)){  
				print "@east_female_nucleotides hello @east_male_nucleotides hello ";
				# this position has two variants in east females population, 
				# and one in east males, so quantify the frequency in other pops. Once concern of doing it this way is that 
				# there could be two homozygous individuals for a different SNP
				# I think this is ok because SNPs often are undercalled.
				$switch=0;
				foreach my $female_nuc (@uniq_east_female_nucleotides){
					foreach my $male_nuc (@uniq_east_male_nucleotides){
						if($female_nuc ne $male_nuc){
							$W_specific = $female_nuc;
							$switch+=1;
						}
					}	
				}
				if(($switch == 1)&&($#uniq_east_female_nucleotides == 1)&&($#uniq_east_male_nucleotides == 0)){ 
					# there was only 1 nucleotide in the females that differed from all the males
					# and all the males had only one nucleotide
					# these are the ones we want to count in the west
					print "$W_specific W ";
					$number_of_alleles=0;
					$number_of_alleles_m=0;
					$number_of_alleles_f=0;
					$count_W_SNPs_in_west_fems=0; # these are counters used for each position
					$count_W_SNPs_in_west_males=0; # these are counters used for each position			  
					foreach my $w_female_nuc (@west_female_nucleotides){ # go through each nucleotide in west females
						if($w_female_nuc eq $W_specific){
							$count_W_SNPs_in_west_fems+=1;
						}
						$number_of_alleles+=1;
						$number_of_alleles_f+=1;	
					}
					foreach my $w_male_nuc (@west_male_nucleotides){ # go through each nucleotide in west males
						if($w_male_nuc eq $W_specific){
							$count_W_SNPs_in_west_males+=1;
						}
						$number_of_alleles+=1;
						$number_of_alleles_m+=1;	
					}
					# now calculate the frequency of that putative W-specific SNP from the east occurs in the west
					if($number_of_alleles > 0){
						$count_W_SNPs_in_west+=1;
						$frequency_of_east_W_SNPs_in_west+=($count_W_SNPs_in_west_fems+$count_W_SNPs_in_west_males)/$number_of_alleles;
					}
					if($number_of_alleles_f > 0){
						$count_W_SNPs_in_west_fems+=1;
						$frequency_of_east_W_SNPs_in_west_fems+=$count_W_SNPs_in_west_fems/$number_of_alleles_f;
					}
					if($number_of_alleles_m > 0){
						$count_W_SNPs_in_west_males+=1;
						$frequency_of_east_W_SNPs_in_west_males+=$count_W_SNPs_in_west_males/$number_of_alleles_m;
					}
					# we will get the average frequency across all W-specific SNPs later by dividing by $count_W_SNPs_in_west
				}
			} # end of check to see if there are only 2 nucleotides in females
		
			# now do the same for the partners of the focal individuals (the males in the east)
			if(($#uniq_east_female_nucleotides == 1)&&($#uniq_east_male_nucleotides == 0)){  
				# this position has two variants in east females, 
				# and one in east males, so quantify the frequency of the putative Z-specific variant 
				# in the west pops. 
				$Z_specific = $uniq_east_male_nucleotides[0];
				print "$Z_specific Z\n";
				$number_of_alleles=0;
				$number_of_alleles_m=0;
				$number_of_alleles_f=0;
				$count_Z_SNPs_in_west_fems=0; # these are counters used for each position
				$count_Z_SNPs_in_west_males=0; # these are counters used for each position		  
				foreach my $w_female_nuc (@west_female_nucleotides){ # go through each nucleotide in west females
					if($w_female_nuc eq $Z_specific){
						$count_Z_SNPs_in_west_fems+=1;
					}
					$number_of_alleles+=1;
					$number_of_alleles_f+=1;	
				}
				foreach my $w_male_nuc (@west_male_nucleotides){ # go through each nucleotide in west males
					if($w_male_nuc eq $Z_specific){
						$count_Z_SNPs_in_west_males+=1;
					}
					$number_of_alleles+=1;
					$number_of_alleles_m+=1;	
				}
				# now calculate the frequency of that putative W-specific SNP from the east occurs in the west
				if($number_of_alleles > 0){
					$count_Z_SNPs_in_west+=1;
					$frequency_of_east_Z_SNPs_in_west+=($count_Z_SNPs_in_west_fems+$count_Z_SNPs_in_west_males)/$number_of_alleles;
				}
				if($number_of_alleles_f > 0){
					$frequency_of_east_Z_SNPs_in_west_fems+=$count_Z_SNPs_in_west_fems/$number_of_alleles_f;
				}
				if($number_of_alleles_m > 0){	
					$frequency_of_east_Z_SNPs_in_west_males+=$count_Z_SNPs_in_west_males/$number_of_alleles_m;
				}
				# we will get the average frequency across all W-specific SNPs later by dividing by $count_Z_SNPs_in_west
			} # end of check to see if there are only 2 nucleotides in females

			
	} # end if to check for header of input file
} # end while	
close DATAINPUT;				

# print out the counts
print "count W_SNPs in west:", $count_W_SNPs_in_west,"\n";
print "count Z_SNPs in west:", $count_Z_SNPs_in_west,"\n\n";
print "Average frequency of east W_SNPs in west:",$frequency_of_east_W_SNPs_in_west/$count_W_SNPs_in_west,"\n";
print "Average frequency of east Z_SNPs in west:",$frequency_of_east_Z_SNPs_in_west/$count_Z_SNPs_in_west,"\n\n";

```
