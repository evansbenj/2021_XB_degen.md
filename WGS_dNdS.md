# Leverage XB genome and XL annotations

One idea is to map the XB data to the XB genome, extract reads from chr8L, and then map these to XL chr8L. This might improve mapping and reduce the amount of cross-paralog mapping. Then the pipeline below could be followed for the female (WZ) and male (ZZ) genomes, most likely using this tree (XL,XB_W,XB_Z) after inferring XB_W and XB_Z by comparing SNPs.


# WGS dNdS

I mapped the WGS data from a female and a male to the XL v9.2 assembly.

I am getting coordinates for each gene on chr8L using this script:
```perl
#!/usr/bin/env perl
use strict;
use warnings;


#  This program reads in gff file and outputs
#  the coordinates of all CDS for each gene, but only for "transcript_01"
#  Other transcripts of each gene are not included; this ensures that each gene is analyzed only once

# to execute type ./Get_coordinates_of_CDS_in_each_gene_outout_individual_beds.pl XLv9.2_xenbase_annotations_chr8L_only.gff
# this inputfile was made using this command (if you use the full gff file, you get coordinates for all genes on all scaffolds and chrs)
# cat XLv9.2_xenbase_annotations.gff | grep 'chr8L' > XLv9.2_xenbase_annotations_chr8L_only.gff

my $inputfile = $ARGV[0];
my $outputfile;

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

my @temp;
my @temp1;
my @Gene_name;
my @Gene_ID;
my @CDS_name;
my @Parent_ID;
my @CDS_Parent_ID;
my $Gene_name;
my $Gene_ID;
my $CDS_name;
my $Parent_ID;
my $CDS_Parent_ID;
my $Transcript_ID;
my %gene_hash;
my $gene_counter=0;
my $exon_counter;
while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	#print $temp[1],"\n";
	if($temp[2] eq 'mRNA'){ # we need to parse the gene name and geneID
		# We are assuming that all CDs are preceded by a "mRNA" annotation
		@temp1=split(';',$temp[8]);
		@Gene_name=split('=',$temp1[0]); 
		$Gene_name=$Gene_name[1]; # for the XL gff file, this begins with "rna" and is the transcript ID for a parent gene
		@Parent_ID=split('=',$temp1[1]); 
		$Parent_ID=$Parent_ID[1]; # for the XL gff file, this begins with "gene" and is the real gene name
		#print $Gene_ID,"\n";
		$gene_counter+=1;
		$exon_counter=1;
		print "hey ",$Parent_ID,"\n";
	}
	elsif($temp[2] eq 'CDS'){ # we need to save the coordinates and phase of this CDS
		# first get the CDS id
		@temp1=split(';',$temp[8]);
		@CDS_name=split('=',$temp1[0]); 
		$CDS_name=$CDS_name[1];# for CDS this is actually a CDS-specific ID
		@CDS_Parent_ID=split('=',$temp1[1]);
		$CDS_Parent_ID=$CDS_Parent_ID[1]; # for CDS this should match the Gene_name of the mRNA row
		print $Gene_name,"\n";
		# only keep CDS that have the same parent ID as the gene name (this is a double check and should work without it)
		if($CDS_Parent_ID eq $Gene_name){
			$gene_hash{$Parent_ID}{$Gene_name}{$exon_counter}[0]=$temp[3]; # start coordinate
			$gene_hash{$Parent_ID}{$Gene_name}{$exon_counter}[1]=$temp[4];	# stop coordinate
			$gene_hash{$Parent_ID}{$Gene_name}{$exon_counter}[2]=$temp[7];	# phase sort of
									#'0' indicates that the first base of the feature is the first base of a codon, 
									#'1' that the second base is the first base of a codon,
									#'2' that the third base is the first base of a codon,
			$gene_hash{$Parent_ID}{$Gene_name}{$exon_counter}[4]=$temp[0]; # chromosome						
			if($temp[6] eq '+'){
				$gene_hash{$Parent_ID}{$Gene_name}{$exon_counter}[3]=1;	# forward orientation

			}
			elsif($temp[6] eq '-'){
				$gene_hash{$Parent_ID}{$Gene_name}{$exon_counter}[3]=-1;	# reverse orientation
			}
			else{
				print "WTF ",$temp[6],"\n";
			}
			$exon_counter+=1;
		}
	}
} # end while	
close DATAINPUT;	
# OK, now all the CDS are in a hash
my $switch=0;
# print a bed file for each transcript (sometimes there will be multiple transcripts for the same gene)
foreach my $Parent_ID (sort keys %gene_hash){
	foreach my $Gene_name (sort keys %{$gene_hash{$Parent_ID}}){
		foreach my $exon (sort {$a <=> $b} keys %{$gene_hash{$Parent_ID}{$Gene_name}}){
			if($switch == 0){ # open a new file for this transcript
				$switch=1;
				if($gene_hash{$Parent_ID}{$Gene_name}{$exon}[3] eq 1){ # assume all exons are in this orientation
					unless (open(OUTFILE, ">".$Parent_ID."_".$Gene_name.".coord"))  {
						print "I can\'t write to $outputfile\n";
						exit;
					}
					print "Creating output file: ".$Parent_ID."_".$Gene_name.".coord\n";
				}
				else{ # assume all exons are in this orientation
					unless (open(OUTFILE, ">".$Parent_ID."_".$Gene_name."_rc.coord"))  {
						print "I can\'t write to $outputfile\n";
						exit;
					}
					print "Creating output file: ".$Parent_ID."_".$Gene_name."_rc.coord\n";
				}	
			}	
			print OUTFILE $gene_hash{$Parent_ID}{$Gene_name}{$exon}[4],"\t";
			print OUTFILE $gene_hash{$Parent_ID}{$Gene_name}{$exon}[0],"\t";
			print OUTFILE $gene_hash{$Parent_ID}{$Gene_name}{$exon}[1],"\n";
		}
		close OUTFILE;
		$switch=0; # reset for next transcript
	}		
}	

```

Now run this command:
```perl
#!/usr/bin/env perl
use strict;
use warnings;


#  This program reads in coordinate files from a directory
# and feeds them into a bash script that extracts sections using bcftools

# to execute type ./Run_bcftools_with_lots_of_inputs.pl path_to_bed_files chr


my $inputfile = $ARGV[0];
	unless (open DATAINPUT, $inputfile) {
		print "Can not find the input file.\n";
		exit;
	}
my $chr = $ARGV[1];
my @files = glob($inputfile.'/*coord');

foreach ( @files ) {
	print $_,"\n";
	system( "./2021_bcftools_extract_sections_from_vcf.sh $_ $chr")
}
```

Which evokes this bash script for each gene:
```sh
#!/bin/sh
#SBATCH --job-name=bcftools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0:10:00
#SBATCH --mem=2gb
#SBATCH --output=bcftools.%J.out
#SBATCH --error=bcftools.%J.err
#SBATCH --account=def-ben

# execute like this: ./2021_bcftools_extract_sections_from_vcf.sh path_and_filename_of_coordinate_file chr
# load these modules before running:
# module load StdEnv/2020 gcc/9.3.0 bcftools/1.11
bcftools view -R ${1} /home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/MandF_WGS_filtered_and_removed_vcfs/al
lsites_${2}.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz -o ${1}.vcf
```

This should generate concatenated vcf files for each gene.

Then I convert it to a tab using vcf2tab.pl:
```perl
#!/usr/bin/env perl
use strict;
use warnings;


# This program makes tab files from each vcf file in a directory

# before executing load modules
# module load nixpkgs/16.09 intel/2018.3 vcftools/0.1.16
# module load StdEnv/2020 perl/5.30.2
# to execute type ./vcf2tab.pl path_to_vcf_files


    
my $inputfile = $ARGV[0];
	unless (open DATAINPUT, $inputfile) {
		print "Can not find the input file.\n";
		exit;
	}

my @files = glob($inputfile.'/*vcf');

foreach ( @files ) {
 #   print $_,"\n";
	system( "vcf-to-tab < $_ > $_.tab")
}
```

Then I convert this to a paml file using tab2paml.pl:
``` perl
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use List::Util qw( min max );


#  This program reads in a tab delimited genotype file generated
#  from vcftools and prints out a paml file 

# Data should be from one chromosome

# on graham, you (may) need to load perl first:
# module load StdEnv/2020 perl/5.30.2

# to execute type tab2paml.pl inputfile.tab output_paml_in 

# From the paml doc:
# In a sequence, T, C, A, G, U, t, c, a, g, u are recognized as nucleotides (for baseml, basemlg and
# codonml), while the standard one-letter codes (A, R, N, D, C, Q, E, G, H, I, L, K, M, F, P, S, T, W,
# Y, V or their lowercase equivalents) are recognized as amino acids. Ambiguity characters
# (undetermined nucleotides or amino acids) are allowed as well. Three special characters ".", "-", and
# "?" are interpreted like this: a dot means the same character as in the first sequence, a dash means an
# alignment gap, and a question mark means an undetermined nucleotide or amino acid.

# if cleandata = 0, both ambiguity characters and alignment gaps are treated
# as ambiguity characters. In the pairwise distance calculation (the lower-diagonal distance matrix in
# the output), cleandata = 1 means “complete deletion”, with all sites involving ambiguity characters
# and alignment gaps removed from all sequences, while cleandata = 0 means “pairwise deletion”,
# with only sites which have missing characters in the pair removed.

# Thus this script will generate one seq for each diploid sequence and use ambiguity characters for heterozygous
# positions, including Ns for deletions or heterozygous deletions

# The analysis will use "cleandata = 0"

# the analysis will include the XL reference seq, a pseudo-Wchr inferred from the female genotypes, 
# and a pseudo-Z chr inferred from the female and male genotypes

# For variable positions in XB, the pseudo-W chr will include female-specific SNPs
# For variable positions in XB, the pseudo-Z chr will include SNPs that are homoz in males but heteroz in females
# For example if the female and male genotypes are C/T and T/T respectively, the W gets the C and the Z gets the T
# if the female and male are heterozygous, the position will be an N

# I'll do the whole chr but it should not be accurate outside of the sex-linked region (>~55Mb)


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

my $XL="";
my $XB_male_SRR6357672_Z="";
my $XB_female_SRR6357673_W="";

my @alleles;
my @temp;
my @lengths;
my $max;
my $a;
my $b;
my $start;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line);
	if($temp[0] ne '#CHROM'){
			# load each of the 4 XL alleles and the XL ref
			for($a=2; $a<=6; $a++){
				$alleles[$a-2]=$temp[$a]; # start index at zero
			}	
			# now see which is the longest
			@lengths=();
			for($a=0; $a<5; $a++){
				push(@lengths,length($alleles[$a]));
			}
			$max = max @lengths;
			#print $max;
			# pad alleles if needed
			my @temp_alleles;
			for($a=0; $a<5; $a++){
				@temp_alleles=split(//,$alleles[$a]);
				if(($#temp_alleles+1) < $max){
					for($b=($#temp_alleles+1); $b<$max; $b++){
						$alleles[$a]=$alleles[$a]."-";
					}	
				}	
			}		
			# now the lengths should all be the same

			# add the data to the chrs
			# first do XL
			$XL=$XL.$alleles[0];
			# now do the XB male, which is the first genotype seq after the ref
			if($alleles[1] eq $alleles[2]){
				$XB_male_SRR6357672_Z=$XB_male_SRR6357672_Z.$alleles[1]; # no polymorphism on the Z
			}
			else{
				for($b=0; $b<$max; $b++){
					$XB_male_SRR6357672_Z=$XB_male_SRR6357672_Z."N"; # heterozygous positions in the male are coded as "N"
				}
			}	
			if($alleles[3] eq $alleles[4]){
				$XB_female_SRR6357673_W=$XB_female_SRR6357673_W.$alleles[3]; # no divergence between W and Z
			}
			else{ # find out which female allele is different from the male allele
				if(($alleles[3] eq $alleles[1])&&($alleles[3] eq $alleles[2])){
					$XB_female_SRR6357673_W=$XB_female_SRR6357673_W.$alleles[4]; # $alleles[4] is the W and $alleles[3] is the Z
				}
				elsif(($alleles[4] eq $alleles[1])&&($alleles[4] eq $alleles[2])){
					$XB_female_SRR6357673_W=$XB_female_SRR6357673_W.$alleles[3]; # $alleles[3] is the W and $alleles[4] is the Z
				}
				else{ # the female alleles are both different from the male alleles or the male is heterozygous
					for($b=0; $b<$max; $b++){
						$XB_female_SRR6357673_W=$XB_female_SRR6357673_W."N";
					}
				}
			}		
			# should be done adding bases to each seq
	} # end if to check for header of input file
} # end while	
close DATAINPUT;				

# print out the lengths of each chr to see if we missed anything
# print "Length of XL :",length($XL),"\n";
# print "Length of Z :",length($XB_male_SRR6357672_Z),"\n";
# print "Length of W :",length($XB_female_SRR6357673_W),"\n";

# reversecomplement if needed
if(substr($inputfile,-17) eq "_rc.coord.vcf.tab"){
	my $XL_revcom=reverse $XL;
	my $XB_male_SRR6357672_Z_revcom=reverse $XB_male_SRR6357672_Z;
	my $XB_female_SRR6357673_W_revcom=reverse $XB_female_SRR6357673_W;
	$XL_revcom =~ tr/ACGTacgt/TGCAtgca/;
	$XB_male_SRR6357672_Z_revcom =~ tr/ACGTacgt/TGCAtgca/;
	$XB_female_SRR6357673_W_revcom =~ tr/ACGTacgt/TGCAtgca/;

	$XL = $XL_revcom;
	$XB_male_SRR6357672_Z = $XB_male_SRR6357672_Z_revcom;
	$XB_female_SRR6357673_W = $XB_female_SRR6357673_W_revcom;
}


# OK print out the fasta file

print OUTFILE2 "3 ",length($XL),"\n";
print OUTFILE2 "XL     ";
print OUTFILE2 $XL,"\n";
print OUTFILE2 "XB_W     ";
print OUTFILE2 $XB_female_SRR6357673_W,"\n";
print OUTFILE2 "XB_Z     ";
print OUTFILE2 $XB_male_SRR6357672_Z,"\n";


close OUTFILE2;
my $n = length($XL)/3;

if( $n != int(length($XL)/3) ){
	print "This infile not multiples of 3: ",$inputfile," ",length($XL),"\n";
}

```
The script above can be evoked for a directory of files using this script:
``` perl
#!/usr/bin/env perl
use strict;
use warnings;


#  This program reads in coordinate files from a directory
# and feeds them into a bash script that extracts sections using bcftools

# before executing load modules
# module load StdEnv/2020 perl/5.30.2
# to execute type ./Make_lots_of_paml_files.pl path_to_tab_files


    
my $inputfile = $ARGV[0];
	unless (open DATAINPUT, $inputfile) {
		print "Can not find the input file.\n";
		exit;
	}

my @files = glob($inputfile.'/*tab');

foreach ( @files ) {
 #   print $_,"\n";
	system( "./tab2paml.pl $_ $_\.paml_in")
}
```
