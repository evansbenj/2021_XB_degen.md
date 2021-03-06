# RADseq data

Are here:
```
/home/ben/projects/rrg-ben/ben/2020_radsex/bin/allsamples
```

# XB RADsex

Previous analyses were done by population on graham here:
```
/home/ben/projects/rrg-ben/ben/2020_radsex/bin
```

I wrote a script to pull out RADsex tags that are present in n males and m females:
```
#!/usr/bin/env perl
use strict;
use warnings;


# This program will extract a fasta file from a RADsex marker file
# that is present in exactly x Males and exactly y Females

# It assumes the sample names begin with "Fem" or "Mal" to indicate the sex

# execute like this:
# ./Gets_RADsex_tags.pl infile n_male_tags n_female_tags

# on computecanada, first load perl: module load StdEnv/2020 perl/5.30.2

my $inputfile1 = $ARGV[0];
my $n_male_tags = $ARGV[1];
my $n_female_tags = $ARGV[2];

unless (open DATAINPUT, $inputfile1) {
	print "Can not find the input file.\n";
	exit;
}

my $outputfile = $inputfile1."_".$n_male_tags."_Males_".$n_female_tags."_Females.fasta"; # the name of the output file is from the commandline
unless (open(OUTFILE, ">$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}

my @sex;
my $counter=0;
my $y;
my $females=0;
my $males=0;
my @temp;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split('\t',$line);
	if(substr($temp[0], 0, 1) ne '#'){ # ignores first line
		if($temp[0] eq 'id'){ 
			for ($y = 0 ; $y < ($#temp+1); $y++ ) {
				if($temp[$y] =~ /^Fem/){
					$sex[$counter]=1; # this is a female
				}
				elsif(substr($temp[$y],0,3) eq "Mal"){
					$sex[$counter]=2; # this is a male
				}
				else{
					$sex[$counter]=0; # this is an index that is not a sample - id or seq
				}
				$counter+=1
			}
			print "hello @sex\n";	
		}
		else{ 
			for ($y = 2 ; $y < ($#sex+1); $y++ ) {
				if($sex[$y] == 1){
					if($temp[$y] > 0){
						$females+=1;
					}
				}
				elsif($sex[$y] == 2){
					if($temp[$y] > 0){
						$males+=1;
					}
				}
				else{
					print "Something weird\n";
				}
			}
			# print $females," ",$males,"\n";
			if(($females==$n_female_tags)&&($males==$n_male_tags)){
				print OUTFILE ">",$temp[0],"\n",$temp[1],"\n";
			}
		}
		$females=0;
		$males=0;
	}	
}		
close DATAINPUT;
close OUTFILE;
```

Let's do this on the Wundyani population first.  The heat plot suggests there are hundreds of tags that are present in zero males and five females but only 5-24 tags that are present in zero females and 5 males. The script identified 228 tags in zero males and five females but only 31 in five males and zero females.

OK let's blast them to the XL genome and summarize the chr and location that they are on.

Load blast
```
module load StdEnv/2020  gcc/9.3.0 blast+/2.11.0
```
Create a blastable db:
```
makeblastdb -in XENLA_9.2_genome.fa -dbtype nucl -out XENLA_9.2_genome.fa_blastable
```
Use the RADsex tags as a blast query:
```
blastn -query wundyanimarkers_table.tsv_0_Males_5_Females.fasta -db /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa_blastable -outfmt 6 -max_target_seqs 1 -out ./wundyanimarkers_table.tsv_0_Males_5_Females.fasta_to_XL_v9.2_genome 

blastn -query wundyanimarkers_table.tsv_5_Males_0_Females.fasta -db /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa_blastable -outfmt 6 -max_target_seqs 1 -out ./wundyanimarkers_table.tsv_5_Males_0_Females.fasta_to_XL_v9.2_genome 

```

Now extract the chr and coord of each hit:
```
cut -f2,9 wundyanimarkers_table.tsv_0_Males_5_Females.fasta_to_XL_v9.2_genome > wundyanimarkers_table.tsv_0_Males_5_Females_hits_to_XL_chr_start
cut -f2,9 wundyanimarkers_table.tsv_5_Males_0_Females.fasta_to_XL_v9.2_genome > wundyanimarkers_table.tsv_5_Males_0_Females_hits_to_XL_chr_start
```

Results are bizarre - female-specific tags don't map to the sex-linked portion of chr8L very much even though we know that this region has very high Fst. This is probably because I only saved the top hit and many tags probably map to repetitive regions.  Instead I'll try saving all hits and only looking at the ones with one match.
```
blastn -query wundyanimarkers_table.tsv_0_Males_5_Females.fasta -db /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa_blastable -outfmt 6 -out ./wundyanimarkers_table.tsv_0_Males_5_Females.fasta_to_XL_v9.2_genome 

blastn -query wundyanimarkers_table.tsv_5_Males_0_Females.fasta -db /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa_blastable -outfmt 6 -out ./wundyanimarkers_table.tsv_5_Males_0_Females.fasta_to_XL_v9.2_genome 

```
Now extract the chr and coord of each hit plus the query ID:
```
cut -f1,2,9 wundyanimarkers_table.tsv_0_Males_5_Females.fasta_to_XL_v9.2_genome > wundyanimarkers_table.tsv_0_Males_5_Females_hits_to_XL_chr_start
cut -f1,2,9 wundyanimarkers_table.tsv_5_Males_0_Females.fasta_to_XL_v9.2_genome > wundyanimarkers_table.tsv_5_Males_0_Females_hits_to_XL_chr_start
```

Now it makes much more sense.  Most of the hits that match only one chr are the sexlinked region of chr8L (RADsextagID top_blasthit_chr top_blast_hit_pos):
```
2687211	chr1L	20661677
2657416	chr2L	108211134
1517301	chr3S	95007302
3553996	chr4L	27940848
1539712	chr5L	140476576
1622854	chr7L	98505191
2116651	chr8L	1104005
2583453	chr8L	2871025
2562926	chr8L	4095162
1619047	chr8L	6047829
3248976	chr8L	6317146
1845941	chr8L	6320724
1592582	chr8L	23170900
2898615	chr8L	30571400
1955550	chr8L	37306450
2121067	chr8L	37957683
3005912	chr8L	43509703
3569224	chr8L	43586650
1623984	chr8L	45286872
1612167	chr8L	46048421
1642858	chr8L	46428193
914951	chr8L	99758497
1529226	chr8S	29279390
1599749	chr8S	53967835
2616460	chr8S	54547471
2122762	chr8S	56807904
1711234	Scaffold158	48922
4347507	Scaffold36	1051525
3132459	Scaffold46	496190
741346	Scaffold65754	192!
```

but only a few random ones for the tags that are present in 5 males but no females:
```
chr7S	17942427	1
chr2L	55974566	1
chr5L	146271619	1

```

Now check the west population to see if these tags are there at all, and if so whether they are female-specific. First get the sequence of each female specific tag on chr8L or chr8S (due to mapping error) below 55Mb:
```
grep -E '2116651|2583453|2562926|1619047|3248976|1845941|1592582|2898615|1955550|2121067|3005912|3569224|1623984|1612167|1642858|1529226|1599749|2616460|2122762' wundyanimarkers_table.tsv | cut -f2
```
Now grep these tags in the western table:
```
head -n2 westkenya_markers_table.tsv > temp
grep -E 'TGCAGGCCAAATCTTGCTTGTGCACTCACTTCCTGTTTCCTGCCTGACAGGAAGTGCAGATTTGGGAGCTCCTGAACCCCATTCATTTCAC|TGCAGGCATGGGGCTACACAGACAGCACTTGCTCAGAAAGTATTTAGGCAACTATTAAGTACATTACATTATCCCCCTGGACTCCCTGTAT|TGCAGGAGGAACACACTGAGCAGCAGAGACAGGGAGGAGCATCCACCTGCTGCAAGCACAGCACTCTCCCTCACTCTCTCTCTCACTTCTT|TGCAGGTACTAGCTGTGCTCTTACCAGTATAACGGAGCTCCGGACCGCTTCAGACACGCTCTGCCTCAACTCCGCCGCGGTTCCCGGTTTG|TGCAGGTTCCCGGAACTCATTTGTTCTGCATTCAGAGCATCGCTCCTTCTCCACTTATGCAAATTACTCTAATGCAGACGAAGCCATAAAC|TGCAGGCTTTATAAATGGATTAGGAAGAAGAAAAGGCTCCATTGGCTTCATTACAGGGAATATATCCCCATTTTTTACTTGAGCTCATGCA|TGCAGGATGTCAAAAGTGATTTAGCTTGGTCTAGTCTATATCAATTAACACAGCTACTATAATAGTGCTGACTGCAGCAGCTGACCTGCTA|TGCAGGGGCTTTACCTGTCAGCTGCATCTGTGCCGAAGAGTCAGTACTGTCAGGGGACTATGGCTTCACCCTGAATGAACGTACAGCCATG|TGCAGGGAGGTTTGGATGACACAGGAATGGTGTTGCACCATTCCACGGAGCAGTGTTGCTTGCACAAACATGATCCACATAAGAGTCATCA|TGCAGGTCTGTCCAATTACGGTGTACAACACACAGCAGAATACGGCACTGCAGCCCCCCAACACGGCAAGGCTCTGGGAGATGACATGGTA|TGCAGGGGAGCAAAAGAGGGTTTCCTGCCTCATGTTTTATGTTGACACTGATTTCTATTTACTCTGCACATTTGAGATGGGAGCTTTGTAG|TGCAGGCCCGGGAAGTCCAATGCCTGGAACCCCTGGTTCTCCCTTTGCACCTGAAACACCTGGAATTCCTGGCTGACCAGGCAGACCTGGT|TGCAGGGGTGCTCAGTAATGTTCTCTAATTTAACCCTTTTCTTATTGCTATGGAAACCTTGAAGCACAAAACAGCTTATTAACTATGCCAG|TGCAGGTAACACTTGTCTGCAGAGGTGCATCAGCAGCTGAGACTGATCAAATCTCCCATCATCCCGTGGGGGAGCATGTATAAAAGGAATG|TGCAGGATTACCCACATCAATCAAGCTCTCTCTAAGGTAAGGCTTTAACATTCAGCCCTTCCGTAGCACAATCCATCAGATCATACCAAGT|TGCAGGCAAGATAAATAAGGATAAATACATGTCATAGGTTTGGGAAATCTCTGTTTTTGTATCTGCATTAGTGCCTAAAAATAACCCCAGT|TGCAGGTTATTTACAGTGTTTCTGTCAAACAGCCCAAGGCCAAAACACAATCTGCTTTGCATGAGTAGTAATGATTTCACGTATTACAGAA|TGCAGGACGCCCTAGGTAAGTCTTCAGCATCTTGCATTTGGTCTTCTGAGGAGTCTGTTTCATTGGCTGTAAGACCTGAAGCTGCTTCTCC|TGCAGGAAAGGCGGTAATAACTATGTTTGTAACATGATGACTTTGCATCCCTCTATAATAACTTCACTTTTAAGTATAATATTTTCAGTAG' westkenya_markers_table.tsv >>temp
```
Now temp should have the markers and counts that are in the west. Out of the 19 female-specific markers from Wundyani, 11 were detected in the west but 2 of these had only 1 tag in one individual, so really 9 out of 19 were detected. Only 1 of 19 was female-specific (3613981) but this one was only detected in two females.  My general interpretation is that this is inconsistent with extinction of the W chromosome and potentially consistent with fixation of the W chromosome.  The other 10 female-specific markers from the east probably are diverged in the west.  This is potentially consistent with extinction of the east Z chromosome in the west population.
```
id	10855	1477937	1745496	1974015	2538998	2609181	2905309	2943208	3345476	3613981
sequence	TGCAGGTACTAGCTGTGCTCTTACCAGTATAACGGAGCTCCGGACCGCTTCAGACACGCTCTGCCTCAACTCCGCCGCGGTTCCCGGTTTG	TGCAGGGAGGTTTGGATGACACAGGAATGGTGTTGCACCATTCCACGGAGCAGTGTTGCTTGCACAAACATGATCCACATAAGAGTCATCA	TGCAGGTCTGTCCAATTACGGTGTACAACACACAGCAGAATACGGCACTGCAGCCCCCCAACACGGCAAGGCTCTGGGAGATGACATGGTA	TGCAGGTAACACTTGTCTGCAGAGGTGCATCAGCAGCTGAGACTGATCAAATCTCCCATCATCCCGTGGGGGAGCATGTATAAAAGGAATG	TGCAGGGGAGCAAAAGAGGGTTTCCTGCCTCATGTTTTATGTTGACACTGATTTCTATTTACTCTGCACATTTGAGATGGGAGCTTTGTAG	TGCAGGATTACCCACATCAATCAAGCTCTCTCTAAGGTAAGGCTTTAACATTCAGCCCTTCCGTAGCACAATCCATCAGATCATACCAAGT	TGCAGGTTCCCGGAACTCATTTGTTCTGCATTCAGAGCATCGCTCCTTCTCCACTTATGCAAATTACTCTAATGCAGACGAAGCCATAAAC	TGCAGGTTATTTACAGTGTTTCTGTCAAACAGCCCAAGGCCAAAACACAATCTGCTTTGCATGAGTAGTAATGATTTCACGTATTACAGAA	TGCAGGCAAGATAAATAAGGATAAATACATGTCATAGGTTTGGGAAATCTCTGTTTTTGTATCTGCATTAGTGCCTAAAAATAACCCCAGT	TGCAGGGGTGCTCAGTAATGTTCTCTAATTTAACCCTTTTCTTATTGCTATGGAAACCTTGAAGCACAAAACAGCTTATTAACTATGCCAG
Fem_Cheplaskei_BJE4459_Xb	0	31	120	63	0	9	4	0	13	0
Fem_Cheplaskei_BJE4461_Xb	0	48	2	56	0	56	37	0	8	0
Fem_Cheplaskei_BJE4470_Xb	0	0	47	59	0	0	51	5	0	5
Fem_Chesuwe_BJE4479_Xb	0	59	9	37	0	0	0	77	0	0
Fem_Chesuwe_BJE4481_Xb	0	0	0	0	0	0	0	0	0	0
Fem_Eldoret_BJE4471_Xb	0	42	13	43	1	0	0	0	0	0
Fem_Eldoret_BJE4472_Xb	0	48	66	0	0	34	19	0	0	0
Fem_Eldoret_BJE4474_Xb	0	0	0	0	0	0	0	0	0	14
Fem_Eldoret_BJE4475_Xb	0	0	0	0	0	0	0	0	0	0
Fem_Eldoret_BJE4476_Xb	0	0	0	0	0	0	0	0	0	0
Fem_Kiminini_BJE4429_Xb	0	0	0	0	0	0	36	0	0	0
Fem_Kiminini_BJE4433_Xb	0	1	47	0	0	0	0	0	31	0
Fem_Lukhome_BJE4441_Xb	0	15	0	0	0	0	0	0	0	0
Fem_Lukhome_BJE4444_Xb	0	0	0	0	0	0	30	0	0	0
Fem_Lukhome_BJE4445_Xb	0	0	91	32	0	0	0	0	25	0
Fem_Lukhome_BJE4446_Xb	0	29	0	53	0	0	0	8	0	0
Mal_Cheplaskei_BJE4460_Xb	0	42	28	0	0	19	0	0	0	0
Mal_Cheplaskei_BJE4462_Xb	0	0	0	0	0	0	0	0	23	0
Mal_Cheplaskei_BJE4465_Xb	0	31	0	4	0	0	0	0	0	0
Mal_Cheplaskei_BJE4469_Xb	0	0	0	0	0	0	0	19	0	0
Mal_Chesuwe_BJE4477_Xb	0	50	35	76	0	57	0	25	33	0
Mal_Chesuwe_BJE4478_Xb	0	0	0	23	0	0	0	0	0	0
Mal_Chesuwe_BJE4480_Xb	0	0	0	0	0	0	0	0	1	0
Mal_Eldoret_BJE4473_Xb	0	57	0	26	0	0	0	0	0	0
Mal_Kiminini_BJE4430_Xb	0	0	42	0	0	0	0	0	0	0
Mal_Kiminini_BJE4431_Xb	0	0	24	14	0	0	0	0	0	0
Mal_Kiminini_BJE4432_Xb	0	0	138	95	0	13	47	0	31	0
Mal_Kisumu_BJE4391_Xb	0	0	0	0	0	0	17	0	0	0
Mal_Lukhome_BJE4442_Xb	0	0	0	0	0	0	0	0	31	0
Mal_Lukhome_BJE4443_Xb	0	52	17	6	0	47	34	0	0	0
Mal_Lukhome_BJE4447_Xb	1	85	136	27	0	101	26	30	0	0

```
I can also search for the other east-female-specific tags in the west by allowing for one mismatch (-1) like this:
```
./agrep/agrep -2 "TGCAGGCCAAATCTTGCTTGTGCACTCACTTCCTGTTTCCTGCCTGACAGGAAGTGCAGATTTGGGAGCTCCTGAACCCCATTCATTTCAC" westkenya_markers_table.tsv > east_female_specific_in_west_one_mismatch
./agrep/agrep -2 "TGCAGGCATGGGGCTACACAGACAGCACTTGCTCAGAAAGTATTTAGGCAACTATTAAGTACATTACATTATCCCCCTGGACTCCCTGTAT" westkenya_markers_table.tsv >> east_female_specific_in_west_one_mismatch
./agrep/agrep -2 "TGCAGGAGGAACACACTGAGCAGCAGAGACAGGGAGGAGCATCCACCTGCTGCAAGCACAGCACTCTCCCTCACTCTCTCTCTCACTTCTT" westkenya_markers_table.tsv >> east_female_specific_in_west_one_mismatch
./agrep/agrep -2 "TGCAGGCTTTATAAATGGATTAGGAAGAAGAAAAGGCTCCATTGGCTTCATTACAGGGAATATATCCCCATTTTTTACTTGAGCTCATGCA" westkenya_markers_table.tsv >> east_female_specific_in_west_one_mismatch
./agrep/agrep -2 "TGCAGGATGTCAAAAGTGATTTAGCTTGGTCTAGTCTATATCAATTAACACAGCTACTATAATAGTGCTGACTGCAGCAGCTGACCTGCTA" westkenya_markers_table.tsv >> east_female_specific_in_west_one_mismatch
./agrep/agrep -2 "TGCAGGGGCTTTACCTGTCAGCTGCATCTGTGCCGAAGAGTCAGTACTGTCAGGGGACTATGGCTTCACCCTGAATGAACGTACAGCCATG" westkenya_markers_table.tsv >> east_female_specific_in_west_one_mismatch
./agrep/agrep -2 "TGCAGGCCCGGGAAGTCCAATGCCTGGAACCCCTGGTTCTCCCTTTGCACCTGAAACACCTGGAATTCCTGGCTGACCAGGCAGACCTGGT" westkenya_markers_table.tsv >> east_female_specific_in_west_one_mismatch
./agrep/agrep -2 "TGCAGGACGCCCTAGGTAAGTCTTCAGCATCTTGCATTTGGTCTTCTGAGGAGTCTGTTTCATTGGCTGTAAGACCTGAAGCTGCTTCTCC" westkenya_markers_table.tsv >> east_female_specific_in_west_one_mismatch
./agrep/agrep -2 "TGCAGGAAAGGCGGTAATAACTATGTTTGTAACATGATGACTTTGCATCCCTCTATAATAACTTCACTTTTAAGTATAATATTTTCAGTAG" westkenya_markers_table.tsv >> east_female_specific_in_west_one_mismatch
```
nomatches with one mismatch

With two mismatches here are the matches:
```
id	1623519	2120724	2732439	1302572	2198206	697584	1792079	1970481
sequence	TGCAGGCTTTATAAATGGATTAGGAAGTAGAAAAGGCTCCATTGGCTTCATTACAGGGAATATATCCCCATTTTTTACTTGAGCTCATGCA	TGCAGGATGTCAAAAGTGACTTAGCTTGGTCTAGTCTATATCAATTAACACAGCTACTATAATAGTGCTGACTGCAGCAGCTGACCTGCTA	TGCAGGGGCTTTACCTGTCAGCTGCATCTGTGCTGAAGAGTCAGTACTGTCAGGGGACTATGGCTTCACCCTGAATGAACGTACAGCCATG	TGCAGGACGCCCTAGGTAAGTCTTCAGCATCTTGCATTTGGTCTTCTGAGGAGTCTGTTTCATTGGCTGTTAGACCTGAAGCTGCTTCTCC	TGCAGGACGCCCTAGGTAAGTCTTCAGCATCTTGCATTTGGTCTTCTGAGGAGTCTGTTTCATTGGCTGTCAGACCTGAAGCTGCTTCTCC	TGCAGGAAAGGCGGTAATAACTATGTTTGAAACATGATGACTTTGCATCCCTCTATAATAACTTCACTTTTAAGTATAATATTTTCAGTAG	TGCAGGAAAGGCGGTAATAACTATGTTTGTAACATGATGACTTTGCATCCCTCTATAATAACTTCACTTTTAGATATAATATTTTCAGTAG	TGCAGGAAAGGCGGTAATAACTATGTTTGTAACATGATGACTTTGCATCCCTCTATAATAACTTCACTTTTAAATATAATATTTTCAGTAG
Fem_Cheplaskei_BJE4459_Xb	0	0	0	0	43	0	0	7
Fem_Cheplaskei_BJE4461_Xb	0	0	0	0	25	0	0	0
Fem_Cheplaskei_BJE4470_Xb	32	65	0	0	0	0	0	49
Fem_Chesuwe_BJE4479_Xb	0	37	4	0	0	0	0	0
Fem_Chesuwe_BJE4481_Xb	0	0	0	0	0	0	0	0
Fem_Eldoret_BJE4471_Xb	0	144	53	0	52	0	0	0
Fem_Eldoret_BJE4472_Xb	21	0	17	0	0	0	0	0
Fem_Eldoret_BJE4474_Xb	0	107	26	0	0	0	0	4
Fem_Eldoret_BJE4475_Xb	0	0	0	0	0	0	0	0
Fem_Eldoret_BJE4476_Xb	0	0	0	0	0	0	0	0
Fem_Kiminini_BJE4429_Xb	0	0	12	0	21	0	0	39
Fem_Kiminini_BJE4433_Xb	0	0	0	1	37	0	0	5
Fem_Lukhome_BJE4441_Xb	0	50	0	0	0	0	0	16
Fem_Lukhome_BJE4444_Xb	0	0	0	0	0	0	0	0
Fem_Lukhome_BJE4445_Xb	39	1	18	0	0	0	0	0
Fem_Lukhome_BJE4446_Xb	0	72	20	0	0	0	0	3
Mal_Cheplaskei_BJE4460_Xb	0	2	0	0	0	0	0	7
Mal_Cheplaskei_BJE4462_Xb	0	0	0	0	0	0	0	10
Mal_Cheplaskei_BJE4465_Xb	0	7	82	0	0	1	0	62
Mal_Cheplaskei_BJE4469_Xb	0	3	0	0	13	0	0	0
Mal_Chesuwe_BJE4477_Xb	0	62	58	0	0	0	0	0
Mal_Chesuwe_BJE4478_Xb	0	0	0	0	0	0	0	21
Mal_Chesuwe_BJE4480_Xb	0	0	47	0	0	0	0	39
Mal_Eldoret_BJE4473_Xb	29	45	43	0	0	0	0	81
Mal_Kiminini_BJE4430_Xb	27	64	0	0	12	0	0	0
Mal_Kiminini_BJE4431_Xb	0	0	37	0	20	0	0	0
Mal_Kiminini_BJE4432_Xb	0	71	26	0	186	0	0	0
Mal_Kisumu_BJE4391_Xb	0	36	0	0	0	0	0	44
Mal_Lukhome_BJE4442_Xb	0	0	0	0	56	0	0	0
Mal_Lukhome_BJE4443_Xb	121	0	45	0	9	0	1	23
Mal_Lukhome_BJE4447_Xb	0	63	72	0	0	0	0	94
```
Three are only in one individual and the others are not female=specific.  4 seqs had no matches with 2 mismatches.3 had 1 matches, 2 had 3 matches.

Also let's do the same thing for the west population to see if any sex-specific tags are clustered in the genome.

First I made a male-specific file like this:
```
ls westkenya_markers_table.tsv_*_Males_0_Females.fasta
	westkenya_markers_table.tsv_10_Males_0_Females.fasta
	westkenya_markers_table.tsv_6_Males_0_Females.fasta
	westkenya_markers_table.tsv_7_Males_0_Females.fasta
	westkenya_markers_table.tsv_8_Males_0_Females.fasta
	westkenya_markers_table.tsv_9_Males_0_Females.fasta
cat westkenya_markers_table.tsv_*_Males_0_Females.fasta > westkenya_markers_table.tsv_6ormore_Males_0_Females.fasta
```
This file has 403 seqs. For comparison, for 6 males 0 females there are 261 seqs but for 0 males and 6 females there are only 47 seqs.  So I think we may have a good signal of a Y chr here, which is expected if the Z chr from the east went extinct.

Blast these putative male-specific seqs against the XL genome:
```
blastn -query westkenya_markers_table.tsv_6ormore_Males_0_Females.fasta -db /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa_blastable -outfmt 6 -out ./westkenya_markers_table.tsv_6ormore_Males_0_Females.fasta_to_XL_v9.2_genome 
```
Now extract the chr and coord of each hit plus the query ID:
```
cut -f1,2,9 westkenya_markers_table.tsv_6ormore_Males_0_Females.fasta_to_XL_v9.2_genome  > westkenya_markers_table.tsv_6ormore_Males_0_Females_hits_to_XL_chr_start
```
Chr7L may be the sex-linked one...although chr8L has several make-specific tags near each other
```
2421308	chr1L	32033673
3197190	chr1L	138394046
2091308	chr1L	172151404
2000168	chr1S	5830731
2104855	chr1S	53258363
3192391	chr1S	132387335
3346637	chr2L	31183530
3992146	chr2L	148209221
2042097	chr2S	53474562
1635440	chr2S	115344896
2076964	chr2S	115344896
2370332	chr2S	115344896
1700204	chr3L	18846092
2264275	chr3L	91708075
1677841	chr3L	117566013
2683802	chr3L	130959633
3786131	chr3L	140187528
2213371	chr3S	55822764
2396532	chr3S	83103495
1588953	chr3S	117760884
2870244	chr4L	13519016
1693147	chr4S	11829414
2912690	chr4S	76569079
1590699	chr5L	107350321
3594888	chr5L	107350321
2591679	chr5S	131341944
3638674	chr6S	26416001
2711466	chr7L	10566181
3953657	chr7L	25048044
1526257	chr7L	36845092
2149042	chr7L	109887398
3930066	chr7L	115767599
1330969	chr7L	124114189
4064444	chr7S	108989
3627445	chr8L	93204786
4002324	chr8L	96663707
3679968	chr8L	108325321
3012046	chr8L	117452680
1726257	chr8S	7034978
3190380	chr8S	47374188
3093469	chr9_10L	16616701
3713635	chr9_10L	106039358
3790404	chr9_10S	35455854
1830252	chr9_10S	38455825
2499804	Scaffold20	564579
3324162	Scaffold28	945006
2529122	Scaffold31	1432730
3006480	Scaffold40682	368
2561908	Scaffold49	334021
2702692	Scaffold49	334021
1618524	Scaffold73	370389
```
Blasts of Scaffold20 to XT suggest it is on chr3. Scaffold28 is probably on chr5. Scaffold31 is ~50Mb in XT on chr3. Scaffold40682 is on chr8 !25Mb in XT. Scaffold49 is on chr5 ~70Mb in XT. Scaffold73 is on chr1 around 215Mb in XT.

# Repeating male-specific tads with Austin's genome

# blast the tads against Austin's genome
```
blastn -query westkenya_markers_table.tsv_6ormore_Males_0_Females.fasta -db ../../2021_Austin_XB_genome/Austin_genome/Xbo.v1.fa_blastable -outfmt 6 -out ./westkenya_markers_table.tsv_6ormore_Males_0_Females.fasta_to_XB
```
Get only columns of interest
```
cut -f1,2,9 westkenya_markers_table.tsv_7ormore_Males_0_Females.fasta_to_XB  > westkenya_markers_table.tsv_7ormore_Males_0_Females_hits_to_XB_chr_start
```
Discard mappings to scaffolds
```
grep 'Chr' westkenya_markers_table.tsv_6ormore_Males_0_Females_hits_to_XB_chr_start > temp
```
Tabulate occurances of first column, which is a RADtag ID
```
cat temp | sort | uniq -c -w7 | sort > tempppp.txt
```
I manually edited tempppp.txt to delete all rows with a first column that was not 1.  Then I plotted this:
```R
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2021_Xborealis_sexchr_degen/RADsex_XB_west')
library (ggplot2)

dat<-read.table("./westkenya_markers_table.tsv_6ormore_Males_0_Females_hits_to_XB_chr_start_chronly_unique_mapping_only",header=F)
# add some headers
names(dat) <- c("one","tag","chr","pos")
# make a column with a random value for plotting
dat$meaningless <- 1
# make a column to color values with low pvalues
dat$color <- "red"


pdf("./West_only_male_specific_AustinGenome.pdf",w=8, h=12.0, version="1.4", bg="transparent")
p<-ggplot(dat, aes(x=pos, y=meaningless, color=color)) + 
    # add points
    geom_point(size=2, alpha = 0.7 ) +
    # add loess line
    # geom_smooth() +
    # label axes
    #labs(x = "XB_W_N_x_dN", y="XB_Z_N_x_dN") +
    # color the stuff the way I want
    scale_colour_manual(name="Genomic region", 
                      values = c("red"="red"),
                      breaks=c("red"),
                      labels=c("male-specific")) +
    facet_wrap(~chr, ncol = 2) +
    # get rid of gray background
    theme_bw() + 
    theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) +
    geom_vline(data = data.frame(xint=140459220,chr="Chr1L"), aes(xintercept = xint), linetype = "dotted")+
    geom_vline(data = data.frame(xint=115397353,chr="Chr1S"), aes(xintercept = xint), linetype = "dotted")+
    geom_vline(data = data.frame(xint=141036519,chr="Chr4L"), aes(xintercept = xint), linetype = "dotted")
p
dev.off()
```
