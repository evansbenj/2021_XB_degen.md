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
