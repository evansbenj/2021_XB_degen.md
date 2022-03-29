# Generate input file for plotting
```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use List::Util qw( min max );



#  This program reads in a tab delimited genotype file generated
#  from vcftools.  It will output a new file with binary genotypes for
#  all biallelic loci: 00, 01, and 11

#  on cedar, you (may) need to load perl first:
#  module load StdEnv/2020 perl/5.30.2

#  to execute type Makes_binary_SNPs_for_haplotype_matrix.pl inputfile.tab output_genotypes 

# Makes_binary_SNPs_for_haplotype_matrix.pl temp.tab temp.out 

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
print "Creating output file: $outputfile\n";

my $counter=2; # this keeps track of the sample IDs
my %genotype_hash;
my $x;
my @temp;
my @lengths;
my @nucleotides;
my @alleles;
my @samples;
my $max;
my $unique;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split(/[\/'\t']+/,$line); # split into individual nucleotides, even within a genotype
	if($temp[0] ne '#CHROM'){
			# load each of the alleles, the first three values are the chr, position, and ref
			@alleles=@temp;	
			$counter=2; # this keeps track of the sample IDs
			# now see which is the longest
			@lengths=();
			for($x=3; $x<=$#alleles; $x++){
				push(@lengths,length($alleles[$x]));
			}
			$max = max @lengths;
			# ignore sites where one or more alleles is longer than 1
			if($max == 1){
				# check if the site is biallelic
				@nucleotides=();
				for($x=3; $x<=$#alleles; $x++){
					if(($alleles[$x] ne ".")&&($alleles[$x] ne "*")){
						push(@nucleotides, $alleles[$x])
					}
				}	
				$unique = uniq @nucleotides;
				if($unique == 2){ # this site is biallelic
					#assign binary genotypes to each individual
					for($x=3; $x<=$#alleles; $x=$x+2){
						$counter+=1;
						if(($alleles[$x] eq $temp[2])&&($alleles[$x+1] eq $temp[2])){ # the genotype is homoz ref
							$genotype_hash{$temp[1]}{$samples[$counter]} = "red"; # the first index is the position and the second index is the sample ID
						}
						elsif(((($alleles[$x] eq $temp[2])&&($alleles[$x+1] ne $temp[2]))||
						(($alleles[$x] ne $temp[2])&&($alleles[$x+1] eq $temp[2])))&&
						($alleles[$x] ne '.')&&($alleles[$x+1] ne '.')&&
						($alleles[$x] ne '*')&&($alleles[$x+1] ne '*'))
						{ # the genotype is heteroz ref
							$genotype_hash{$temp[1]}{$samples[$counter]} = "purple"; # the first index is the position and the second index is the sample ID
						}
						elsif((($alleles[$x] ne $temp[2])&&($alleles[$x+1] ne $temp[2]))&&
						($alleles[$x] ne '.')&&($alleles[$x+1] ne '.')&&
						($alleles[$x] ne '*')&&($alleles[$x+1] ne '*'))
						{ # the genotype is heteroz ref
							$genotype_hash{$temp[1]}{$samples[$counter]} = "blue"; # the first index is the position and the second index is the sample ID
						}
						else{ # the genotype is something else - missing 
							$genotype_hash{$temp[1]}{$samples[$counter]} = "white";
						}
					}	
				}	
			}
	} # end if to check for header of input file
	else{
		# record sample names
		@samples = @temp;
		print $#samples,"\n";
	}
} # end while	
close DATAINPUT;				

# print out the counts
print OUTFILE "pos\t";
for($x=3; $x<=$#samples; $x++){
	print OUTFILE $samples[$x],"\t";
}
print OUTFILE "\n";	
# keys of the outer hash
foreach my $position (sort { $a <=> $b } keys %genotype_hash) {
    # keys of the inner hash
	print OUTFILE $position,"\t";
    foreach my $sample (sort keys %{$genotype_hash{$position}}) {
        print OUTFILE $genotype_hash{$position}{$sample},"\t";
    }
	print OUTFILE "\n";
}
close OUTFILE;

```

# Plot
```R
library(dplyr)
library(ggplot2)
library("reshape") 

setwd('/Users/Shared/Previously Relocated Items/Security/projects/submitted/2021_Xborealis_sexchr_degen/RADseq_counts_in_west_of_SNPs_in_east')
dat<-read.table("./temp.out",header=TRUE)

# subset to include only Wundyani individuals
wundyani <- dat[,c(1,26:30,51:55)]
tail(wundyani)
wundyani_SL <- filter(wundyani, pos < 105733990) # length 105733990
                                                # or for SL region 23000000
dim(wundyani_SL)
wundyani_SL$pos <- factor(wundyani_SL$pos)
wundyani_SL_sorted <- wundyani_SL[order(as.numeric(as.character(wundyani_SL$pos))), ]

rows_to_remove <- NULL
# Remove rows that are invariant
for(i in 1:nrow(wundyani_SL_sorted)){
    vector <- NULL
    for(j in 2:ncol(wundyani_SL_sorted)){
       vector<- append(vector,as.character(wundyani_SL_sorted[i,j]))
    } 
    if(length(unique(vector)) == 1){
        # get rid of this row because it is invariant
        rows_to_remove <- append(rows_to_remove,i)
    }
    if((length(unique(vector)) == 2)&&('white' %in% vector)){
        # get rid of this row because it is also invariant
        # because of missing data
        rows_to_remove <- append(rows_to_remove,i)
    }
}
dim(wundyani_SL_sorted)
wundyani_SL_sorted <- wundyani_SL_sorted[-rows_to_remove, ]
dim(wundyani_SL_sorted)
#wundyani_SL <- filter(wundyani, pos < 23000000)
names(wundyani_SL_sorted)[2] <- "Fem BJE4515"
names(wundyani_SL_sorted)[3] <- "Fem BJE4516"
names(wundyani_SL_sorted)[4] <- "Fem BJE4534"
names(wundyani_SL_sorted)[5] <- "Fem BJE4535"
names(wundyani_SL_sorted)[6] <- "Fem BJE4541"
names(wundyani_SL_sorted)[7] <- "Mal BJE4536"
names(wundyani_SL_sorted)[8] <- "Mal BJE4537"
names(wundyani_SL_sorted)[9] <- "Mal BJE4538"
names(wundyani_SL_sorted)[10] <- "Mal BJE4539"
names(wundyani_SL_sorted)[11] <- "Mal BJE4540"

#heatmap(as.matrix(wundyani), Rowv = NA, Colv = NA)
#wundyani_SL_sorted$pos <- as.character(wundyani_SL_sorted$pos)


 
wundyani_SL_sorted_melt <- melt(wundyani_SL_sorted,id=c("pos"))    # Reorder data
head(wundyani_SL_sorted_melt)

# sort the taxa
wundyani_SL_sorted_melt$variable <- factor(wundyani_SL_sorted_melt$variable, 
                                 ordered=TRUE, 
                                 levels = c("Fem BJE4516","Fem BJE4541","Fem BJE4535", 
                                            "Mal BJE4536","Mal BJE4539", 
                                            "Fem BJE4515", "Fem BJE4534","Mal BJE4537",
                                            "Mal BJE4538","Mal BJE4540"))

# sort the positions
wundyani_SL_sorted_melt$pos <- factor(wundyani_SL_sorted_melt$pos, 
                            ordered=TRUE)
wundyani_SL_sorted_melt <- wundyani_SL_sorted_melt[order(as.numeric(as.character(wundyani_SL_sorted_melt$pos))), ]
#head(wundyani)
#wundyani_melt$value <- as.character(wundyani_melt$value)
dim(wundyani_SL_sorted_melt)
names(wundyani_SL_sorted_melt)[1] <- "Position"
names(wundyani_SL_sorted_melt)[2] <- "Sample"
ggp <- ggplot(wundyani_SL_sorted_melt, aes(Position, Sample)) +  
    geom_tile(aes(fill = value))+ 
    scale_fill_manual(values=c(red="light blue", purple = "dark gray", blue="black", 
                 white = "white"),
                 labels = c(red="homoz ref", purple = "heteroz", 
                            blue="homoz alt", white = "missing")) +
    theme_classic() +
    #theme(axis.text.x=element_text(angle=45,hjust=1)) +
    theme(axis.text.x=element_blank(),axis.ticks.x = element_blank()) +
    geom_rect(data=wundyani_SL_sorted_melt, 
              mapping=aes(xmin=2, xmax=195, 
                          ymin=7.5, ymax=10.5), 
              color="red", alpha=0, size=1) +
    geom_rect(data=wundyani_SL_sorted_melt, 
              mapping=aes(xmin=2, xmax=195, 
                          ymin=5.5, ymax=7.4), 
              color="blue", alpha=0, size=1) #+
#    geom_rect(data=wundyani_SL_sorted_melt, 
#              mapping=aes(xmin=2, xmax=195, 
#                          ymin=0.5, ymax=5.4), 
#              color="black", alpha=0, size=1)

pdf("./Chr7S_haplotypeMatrix.pdf",w=8, h=2.0, version="1.4", bg="transparent")
    ggp   
dev.off()


```
