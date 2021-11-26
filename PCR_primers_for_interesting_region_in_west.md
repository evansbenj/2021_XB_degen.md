# Interesting region in west  
The WGS data identified an interesting region in the west.  I am going to design primers that target CDS that have het positions in all females but no males.

To do this I need to get coordinates from XLv10 CDS, extract their seqs, and blast them against Austin's genome to get the XB sequences.

I got the CDS coordinates from the XLv10 gff file and got the fast seqs using bedtools:

Working on graham in this directory:
```
/home/ben/projects/rrg-ben/ben/2021_XL_v10_refgenome
```
Get fasta seq:
```
module load StdEnv/2020 bedtools/2.30.0
bedtools getfasta -fi XENLA_10.1_genome.fa -bed CDS_in_Chr8L_XXXX_Mb.bed -fo XL_v10_Chr8L_XXX_Mb_CDSonly.fasta
```

Now blast them against Austin genome:
```
module load StdEnv/2020  gcc/9.3.0 blast+/2.11.0
blastn -query XL_v10_Chr8L_XXX_Mb_CDSonly.fasta -db Xbo.v1.fa_blastable -outfmt 6 -out ./XL_v10_Chr8L_XXX_Mb_CDSonly.fasta_to_Austin_genome 

```
