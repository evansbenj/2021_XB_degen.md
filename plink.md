# Plink

Associations between SNPs and a phenotype (such as phenotypic sex) are easily calculated using plink.

```
module load nixpkgs/16.09 plink/1.9b_5.2-x86_64

plink --vcf temp.vcf.gz --recode --const-fid 0 --chr-set 36 no-y no-xy no-mt --out myplink
plink --file myplink --pheno sex_phenotype --assoc --allow-no-sex
```
where the "sex_phenotype" file is a tab-delimited file that looks like this:
```
0	sample1	1
0	sample2	2
0	sample3	1
0	sample4	2
0	sample5	2
0	sample6	1
```
The first column is the family ID (just zeros here).  The second column is the sample name - this is the same as in the vcf file.  The third column is the phenotype - should use 1 and 2, NOT 1 and 0, because 0 *might* be interpreted as a missing phenotye.

Also this flag "--const-fid 0" sets the family id to zero and tells plink to use the vcf sample name as the sample ID irrespective of whether there is an underscore in the name.

This flag "--chr-set 36" allows extra chrs.  They will be numbers in the order they are encountered in the vcf file (I think - this will need to be confirmed...)

# Subsetting vcf

Working in this directory:
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome
```
concatenate:
```
module load StdEnv/2020 vcftools/0.1.16
vcf-concat DB_chr1L_out.vcf DB_chr1S_out.vcf DB_chr2L_out.vcf DB_chr2S_out.vcf DB_chr3L_out.vcf DB_chr3S_out.vcf DB_chr4L_out.vcf DB_chr4S_out.vcf DB_chr5L_out.vcf DB_chr5S_out.vcf DB_chr6L_out.vcf DB_chr6S_out.vcf DB_chr7L_out.vcf DB_chr7S_out.vcf DB_chr8L_out.vcf DB_chr8S_out.vcf DB_chr9_10L_out.vcf DB_chr9_10S_out.vcf | gzip -c > XB_unfiltered_allchrs.vcf.gz
```
Or filter to include only chrs for XB_Au
```
module load StdEnv/2020  gcc/9.3.0 bcftools/1.11

bcftools view mpileup_raw_wildBorealis_AuXXXGenome.vcf.gz --regions Chr1L,Chr1S,Chr2L,Chr2S,Chr3L,Chr3S,Chr4L,Chr4S,Chr5L,Chr5S,Chr6L,Chr6S,Chr7L,Chr7S,Chr8L,Chr8S,Chr9_10L,Chr9_10S -O z -o mpileup_raw_wildBorealis_AuXXXGenome_Chrs_only.vcf.gz
```

Subset populations
```
vcf-subset -c Fem_Cheplaskei_BJE4459_Xb.fq.gz,Fem_Cheplaskei_BJE4461_Xb.fq.gz,Fem_Cheplaskei_BJE4470_Xb.fq.gz,Fem_Chesuwe_BJE4479_Xb.fq.gz,Fem_Chesuwe_BJE4481_Xb.fq.gz,Fem_Eldoret_BJE4471_Xb.fq.gz,Fem_Eldoret_BJE4472_Xb.fq.gz,Fem_Eldoret_BJE4474_Xb.fq.gz,Fem_Eldoret_BJE4475_Xb.fq.gz,Fem_Eldoret_BJE4476_Xb.fq.gz,Fem_Kiminini_BJE4429_Xb.fq.gz,Fem_Kiminini_BJE4433_Xb.fq.gz,Fem_Lukhome_BJE4441_Xb.fq.gz,Fem_Lukhome_BJE4444_Xb.fq.gz,Fem_Lukhome_BJE4445_Xb.fq.gz,Fem_Lukhome_BJE4446_Xb.fq.gz,Mal_Cheplaskei_BJE4460_Xb.fq.gz,Mal_Cheplaskei_BJE4462_Xb.fq.gz,Mal_Cheplaskei_BJE4465_Xb.fq.gz,Mal_Cheplaskei_BJE4469_Xb.fq.gz,Mal_Chesuwe_BJE4477_Xb.fq.gz,Mal_Chesuwe_BJE4478_Xb.fq.gz,Mal_Chesuwe_BJE4480_Xb.fq.gz,Mal_Eldoret_BJE4473_Xb.fq.gz,Mal_Kiminini_BJE4430_Xb.fq.gz,Mal_Kiminini_BJE4431_Xb.fq.gz,Mal_Kiminini_BJE4432_Xb.fq.gz,Mal_Kisumu_BJE4391_Xb.fq.gz,Mal_Lukhome_BJE4442_Xb.fq.gz,Mal_Lukhome_BJE4443_Xb.fq.gz,Mal_Lukhome_BJE4447_Xb.fq.gz XB_unfiltered_allchrs.vcf.gz | bgzip -c > XB_westonly_unfiltered_allchrs.vcf.gz

vcf-subset -c Fem_NMK_BJE4562_Xb.fq.gz,Fem_NMK_BJE4564_Xb.fq.gz,Fem_NMK_BJE4567_Xb.fq.gz,Fem_NMK_BJE4568_Xb.fq.gz,Fem_Nakuru_BJE4363_Xb.fq.gz,Fem_Nakuru_BJE4364_Xb.fq.gz,Fem_Nakuru_BJE4367_Xb.fq.gz,Fem_Nakuru_BJE4368_Xb.fq.gz,Mal_NMK_BJE4563_Xb.fq.gz,Mal_Nakuru_BJE4365_Xb.fq.gz,Mal_Nakuru_BJE4366_Xb.fq.gz,Mal_Nakuru_BJE4374_Xb.fq.gz,Mal_Nakuru_BJE4377_Xb.fq.gz XB_unfiltered_allchrs.vcf.gz | bgzip -c > XB_centralonly_unfiltered_allchrs.vcf.gz

vcf-subset -c Fem_Wundanyi_BJE4515_Xb.fq.gz,Fem_Wundanyi_BJE4516_Xb.fq.gz,Fem_Wundanyi_BJE4534_Xb.fq.gz,Fem_Wundanyi_BJE4535_Xb.fq.gz,Fem_Wundanyi_BJE4541_Xb.fq.gz,Mal_Wundanyi_BJE4536_Xb.fq.gz,Mal_Wundanyi_BJE4537_Xb.fq.gz,Mal_Wundanyi_BJE4538_Xb.fq.gz,Mal_Wundanyi_BJE4539_Xb.fq.gz,Mal_Wundanyi_BJE4540_Xb.fq.gz XB_unfiltered_allchrs.vcf.gz | bgzip -c > XB_east_Wudyani_only_unfiltered_allchrs.vcf.gz

vcf-subset -c Fem_Nakuru_BJE4363_Xb.fq.gz,Fem_Nakuru_BJE4364_Xb.fq.gz,Fem_Nakuru_BJE4367_Xb.fq.gz,Fem_Nakuru_BJE4368_Xb.fq.gz,Mal_Nakuru_BJE4365_Xb.fq.gz,Mal_Nakuru_BJE4366_Xb.fq.gz,Mal_Nakuru_BJE4374_Xb.fq.gz,Mal_Nakuru_BJE4377_Xb.fq.gz XB_unfiltered_allchrs.vcf.gz | bgzip -c > XB_Njoroonly_unfiltered_allchrs.vcf.gz

vcf-subset -c Fem_Cheplaskei_BJE4459_Xb.fq.gz,Fem_Cheplaskei_BJE4461_Xb.fq.gz,Fem_Cheplaskei_BJE4470_Xb.fq.gz,Fem_Chesuwe_BJE4479_Xb.fq.gz,Fem_Chesuwe_BJE4481_Xb.fq.gz,Fem_Eldoret_BJE4471_Xb.fq.gz,Fem_Eldoret_BJE4472_Xb.fq.gz,Fem_Eldoret_BJE4474_Xb.fq.gz,Fem_Eldoret_BJE4475_Xb.fq.gz,Fem_Eldoret_BJE4476_Xb.fq.gz,Fem_Kiminini_BJE4429_Xb.fq.gz,Fem_Kiminini_BJE4433_Xb.fq.gz,Fem_Lukhome_BJE4441_Xb.fq.gz,Fem_Lukhome_BJE4444_Xb.fq.gz,Fem_Lukhome_BJE4445_Xb.fq.gz,Fem_Lukhome_BJE4446_Xb.fq.gz,Fem_Nakuru_BJE4363_Xb.fq.gz,Fem_Nakuru_BJE4364_Xb.fq.gz,Fem_Nakuru_BJE4367_Xb.fq.gz,Fem_Nakuru_BJE4368_Xb.fq.gz,Mal_Cheplaskei_BJE4460_Xb.fq.gz,Mal_Cheplaskei_BJE4462_Xb.fq.gz,Mal_Cheplaskei_BJE4465_Xb.fq.gz,Mal_Cheplaskei_BJE4469_Xb.fq.gz,Mal_Chesuwe_BJE4477_Xb.fq.gz,Mal_Chesuwe_BJE4478_Xb.fq.gz,Mal_Chesuwe_BJE4480_Xb.fq.gz,Mal_Eldoret_BJE4473_Xb.fq.gz,Mal_Kiminini_BJE4430_Xb.fq.gz,Mal_Kiminini_BJE4431_Xb.fq.gz,Mal_Kiminini_BJE4432_Xb.fq.gz,Mal_Kisumu_BJE4391_Xb.fq.gz,Mal_Lukhome_BJE4442_Xb.fq.gz,Mal_Lukhome_BJE4443_Xb.fq.gz,Mal_Lukhome_BJE4447_Xb.fq.gz,Mal_Nakuru_BJE4365_Xb.fq.gz,Mal_Nakuru_BJE4366_Xb.fq.gz,Mal_Nakuru_BJE4374_Xb.fq.gz,Mal_Nakuru_BJE4377_Xb.fq.gz XB_unfiltered_allchrs.vcf.gz | bgzip -c > XB_west_and_Njoro_only_unfiltered_allchrs.vcf.gz

```
Or use bcftools (much faster):
```
module load StdEnv/2020  gcc/9.3.0 bcftools/1.13

bcftools view -s Fem_Cheplaskei_BJE4459_Xb_trim_sorted.bam,Fem_Cheplaskei_BJE4461_Xb_trim_sorted.bam,Fem_Cheplaskei_BJE4470_Xb_trim_sorted.bam,Fem_Chesuwe_BJE4479_Xb_trim_sorted.bam,Fem_Chesuwe_BJE4481_Xb_trim_sorted.bam,Fem_Eldoret_BJE4471_Xb_trim_sorted.bam,Fem_Eldoret_BJE4472_Xb_trim_sorted.bam,Fem_Eldoret_BJE4474_Xb_trim_sorted.bam,Fem_Eldoret_BJE4475_Xb_trim_sorted.bam,Fem_Eldoret_BJE4476_Xb_trim_sorted.bam,Fem_Kiminini_BJE4429_Xb_trim_sorted.bam,Fem_Kiminini_BJE4433_Xb_trim_sorted.bam,Fem_Lukhome_BJE4441_Xb_trim_sorted.bam,Fem_Lukhome_BJE4444_Xb_trim_sorted.bam,Fem_Lukhome_BJE4445_Xb_trim_sorted.bam,Fem_Lukhome_BJE4446_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4460_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4462_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4465_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4469_Xb_trim_sorted.bam,Mal_Chesuwe_BJE4477_Xb_trim_sorted.bam,Mal_Chesuwe_BJE4478_Xb_trim_sorted.bam,Mal_Chesuwe_BJE4480_Xb_trim_sorted.bam,Mal_Eldoret_BJE4473_Xb_trim_sorted.bam,Mal_Kiminini_BJE4430_Xb_trim_sorted.bam,Mal_Kiminini_BJE4431_Xb_trim_sorted.bam,Mal_Kiminini_BJE4432_Xb_trim_sorted.bam,Mal_Kisumu_BJE4391_Xb_trim_sorted.bam,Mal_Lukhome_BJE4442_Xb_trim_sorted.bam,Mal_Lukhome_BJE4443_Xb_trim_sorted.bam,Mal_Lukhome_BJE4447_Xb_trim_sorted.bam XBwild_unfiltered_allchrs.vcf.gz -O z > west_only_XBwild_unfiltered_allchrs.vcf.gz

bcftools view -s Fem_NMK_BJE4562_Xb_trim_sorted.bam,Fem_NMK_BJE4564_Xb_trim_sorted.bam,Fem_NMK_BJE4567_Xb_trim_sorted.bam,Fem_NMK_BJE4568_Xb_trim_sorted.bam,Mal_NMK_BJE4563_Xb_trim_sorted.bam XBwild_unfiltered_allchrs.vcf.gz -O z > NMK_only_XBwild_unfiltered_allchrs.vcf.gz


bcftools view -s Fem_Cheplaskei_BJE4459_Xb_trim_sorted.bam,Fem_Cheplaskei_BJE4461_Xb_trim_sorted.bam,Fem_Cheplaskei_BJE4470_Xb_trim_sorted.bam,Fem_Chesuwe_BJE4479_Xb_trim_sorted.bam,Fem_Chesuwe_BJE4481_Xb_trim_sorted.bam,Fem_Eldoret_BJE4471_Xb_trim_sorted.bam,Fem_Eldoret_BJE4472_Xb_trim_sorted.bam,Fem_Eldoret_BJE4474_Xb_trim_sorted.bam,Fem_Eldoret_BJE4475_Xb_trim_sorted.bam,Fem_Eldoret_BJE4476_Xb_trim_sorted.bam,Fem_Kiminini_BJE4429_Xb_trim_sorted.bam,Fem_Kiminini_BJE4433_Xb_trim_sorted.bam,Fem_Lukhome_BJE4441_Xb_trim_sorted.bam,Fem_Lukhome_BJE4444_Xb_trim_sorted.bam,Fem_Lukhome_BJE4445_Xb_trim_sorted.bam,Fem_Lukhome_BJE4446_Xb_trim_sorted.bam,Fem_Nakuru_BJE4363_Xb_trim_sorted.bam,Fem_Nakuru_BJE4364_Xb_trim_sorted.bam,Fem_Nakuru_BJE4367_Xb_trim_sorted.bam,Fem_Nakuru_BJE4368_Xb_trim_sorted.bam,Fem_NMK_BJE4562_Xb_trim_sorted.bam,Fem_NMK_BJE4564_Xb_trim_sorted.bam,Fem_NMK_BJE4567_Xb_trim_sorted.bam,Fem_NMK_BJE4568_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4460_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4462_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4465_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4469_Xb_trim_sorted.bam,Mal_Chesuwe_BJE4477_Xb_trim_sorted.bam,Mal_Chesuwe_BJE4478_Xb_trim_sorted.bam,Mal_Chesuwe_BJE4480_Xb_trim_sorted.bam,Mal_Eldoret_BJE4473_Xb_trim_sorted.bam,Mal_Kiminini_BJE4430_Xb_trim_sorted.bam,Mal_Kiminini_BJE4431_Xb_trim_sorted.bam,Mal_Kiminini_BJE4432_Xb_trim_sorted.bam,Mal_Kisumu_BJE4391_Xb_trim_sorted.bam,Mal_Lukhome_BJE4442_Xb_trim_sorted.bam,Mal_Lukhome_BJE4443_Xb_trim_sorted.bam,Mal_Lukhome_BJE4447_Xb_trim_sorted.bam,Mal_Nakuru_BJE4365_Xb_trim_sorted.bam,Mal_Nakuru_BJE4366_Xb_trim_sorted.bam,Mal_Nakuru_BJE4374_Xb_trim_sorted.bam,Mal_Nakuru_BJE4377_Xb_trim_sorted.bam,Mal_NMK_BJE4563_Xb_trim_sorted.bam XBwild_unfiltered_allchrs.vcf.gz -O z > westandNMKandNakuru_XBwild_unfiltered_allchrs.vcf.gz



bcftools view -s Fem_Nakuru_BJE4363_Xb_trim_sorted.bam,Fem_Nakuru_BJE4364_Xb_trim_sorted.bam,Fem_Nakuru_BJE4367_Xb_trim_sorted.bam,Fem_Nakuru_BJE4368_Xb_trim_sorted.bam,Mal_Nakuru_BJE4365_Xb_trim_sorted.bam,Mal_Nakuru_BJE4366_Xb_trim_sorted.bam,Mal_Nakuru_BJE4374_Xb_trim_sorted.bam,Mal_Nakuru_BJE4377_Xb_trim_sorted.bam XBwild_unfiltered_allchrs.vcf.gz -O z > Nakuru_only_XBwild_unfiltered_allchrs.vcf.gz

bcftools view -s Fem_Cheplaskei_BJE4459_Xb_trim_sorted.bam,Fem_Cheplaskei_BJE4461_Xb_trim_sorted.bam,Fem_Cheplaskei_BJE4470_Xb_trim_sorted.bam,Fem_Chesuwe_BJE4479_Xb_trim_sorted.bam,Fem_Chesuwe_BJE4481_Xb_trim_sorted.bam,Fem_Eldoret_BJE4471_Xb_trim_sorted.bam,Fem_Eldoret_BJE4472_Xb_trim_sorted.bam,Fem_Eldoret_BJE4474_Xb_trim_sorted.bam,Fem_Eldoret_BJE4475_Xb_trim_sorted.bam,Fem_Eldoret_BJE4476_Xb_trim_sorted.bam,Fem_Kiminini_BJE4429_Xb_trim_sorted.bam,Fem_Kiminini_BJE4433_Xb_trim_sorted.bam,Fem_Lukhome_BJE4441_Xb_trim_sorted.bam,Fem_Lukhome_BJE4444_Xb_trim_sorted.bam,Fem_Lukhome_BJE4445_Xb_trim_sorted.bam,Fem_Lukhome_BJE4446_Xb_trim_sorted.bam,Fem_Nakuru_BJE4363_Xb_trim_sorted.bam,Fem_Nakuru_BJE4364_Xb_trim_sorted.bam,Fem_Nakuru_BJE4367_Xb_trim_sorted.bam,Fem_Nakuru_BJE4368_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4460_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4462_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4465_Xb_trim_sorted.bam,Mal_Cheplaskei_BJE4469_Xb_trim_sorted.bam,Mal_Chesuwe_BJE4477_Xb_trim_sorted.bam,Mal_Chesuwe_BJE4478_Xb_trim_sorted.bam,Mal_Chesuwe_BJE4480_Xb_trim_sorted.bam,Mal_Eldoret_BJE4473_Xb_trim_sorted.bam,Mal_Kiminini_BJE4430_Xb_trim_sorted.bam,Mal_Kiminini_BJE4431_Xb_trim_sorted.bam,Mal_Kiminini_BJE4432_Xb_trim_sorted.bam,Mal_Kisumu_BJE4391_Xb_trim_sorted.bam,Mal_Lukhome_BJE4442_Xb_trim_sorted.bam,Mal_Lukhome_BJE4443_Xb_trim_sorted.bam,Mal_Lukhome_BJE4447_Xb_trim_sorted.bam,Mal_Nakuru_BJE4365_Xb_trim_sorted.bam,Mal_Nakuru_BJE4366_Xb_trim_sorted.bam,Mal_Nakuru_BJE4374_Xb_trim_sorted.bam,Mal_Nakuru_BJE4377_Xb_trim_sorted.bam XBwild_unfiltered_allchrs.vcf.gz -O z > westandNakuru_XBwild_unfiltered_allchrs.vcf.gz

bcftools view -s Fem_NMK_BJE4562_Xb_trim_sorted.bam,Fem_NMK_BJE4564_Xb_trim_sorted.bam,Fem_NMK_BJE4567_Xb_trim_sorted.bam,Fem_NMK_BJE4568_Xb_trim_sorted.bam,Fem_Nakuru_BJE4363_Xb_trim_sorted.bam,Fem_Nakuru_BJE4364_Xb_trim_sorted.bam,Fem_Nakuru_BJE4367_Xb_trim_sorted.bam,Fem_Nakuru_BJE4368_Xb_trim_sorted.bam,Mal_Nakuru_BJE4365_Xb_trim_sorted.bam,Mal_Nakuru_BJE4366_Xb_trim_sorted.bam,Mal_Nakuru_BJE4374_Xb_trim_sorted.bam,Mal_Nakuru_BJE4377_Xb_trim_sorted.bam,Mal_NMK_BJE4563_Xb_trim_sorted.bam XBwild_unfiltered_allchrs.vcf.gz -O z > NMKandNakuru_only_XBwild_unfiltered_allchrs.vcf.gz

bcftools view -s Fem_Wundanyi_BJE4515_Xb_trim_sorted.bam,Fem_Wundanyi_BJE4516_Xb_trim_sorted.bam,Fem_Wundanyi_BJE4534_Xb_trim_sorted.bam,Fem_Wundanyi_BJE4535_Xb_trim_sorted.bam,Fem_Wundanyi_BJE4541_Xb_trim_sorted.bam,Mal_Wundanyi_BJE4536_Xb_trim_sorted.bam,Mal_Wundanyi_BJE4537_Xb_trim_sorted.bam,Mal_Wundanyi_BJE4538_Xb_trim_sorted.bam,Mal_Wundanyi_BJE4539_Xb_trim_sorted.bam,Mal_Wundanyi_BJE4540_Xb_trim_sorted.bam XBwild_unfiltered_allchrs.vcf.gz -O z > east_only_XBwild_unfiltered_allchrs.vcf.gz
```


# Plink
```
module load nixpkgs/16.09 plink/1.9b_5.2-x86_64

plink --vcf XB_westonly_unfiltered_allchrs.vcf.gz --recode --const-fid 0 --chr-set 36 no-y no-xy no-mt --allow-extra-chr --out XB_westonly_unfiltered_allchrs.vcf.gz_myplink
plink --file XB_westonly_unfiltered_allchrs.vcf.gz_myplink --pheno west_sample_sex --assoc --allow-no-sex --allow-extra-chr 
mv plink.assoc XB_westonly_unfiltered_allchrs.vcf.gz_plink.assoc

plink --vcf XB_Njoroonly_unfiltered_allchrs.vcf.gz --recode --const-fid 0 --chr-set 36 no-y no-xy no-mt --allow-extra-chr --out XB_Njoro_only_unfiltered_allchrs.vcf.gz_myplink
plink --file XB_Njoro_only_unfiltered_allchrs.vcf.gz_myplink --pheno Njoro_sample_sex --assoc --allow-no-sex --allow-extra-chr 
mv plink.assoc XB_Njoro_only_unfiltered_allchrs.vcf.gz_plink.assoc

plink --vcf XB_west_and_Njoro_only_unfiltered_allchrs.vcf.gz --recode --const-fid 0 --chr-set 36 no-y no-xy no-mt --allow-extra-chr --out XB_west_and_Njoro_only_unfiltered_allchrs.vcf.gz_myplink
plink --file XB_west_and_Njoro_only_unfiltered_allchrs.vcf.gz_myplink --pheno west_and_Njoro_sample_sex --assoc --allow-no-sex --allow-extra-chr 
mv plink.assoc XB_west_and_Njoro_only_unfiltered_allchrs.vcf.gz_plink.assoc
```

# Checking for specificity of SL regions
For most plots (except from the east), there is a fair amount of background signal. This makes it hard to know what actually might be sex-linked.  To try to tease this apart, I am going to extract genomic windows that surround the sex-linked SNPs, and blast these against the genome to identify ones with unique matches.  Then I can reassess what to believe.

First make a bedfile from the plink output that had the SNPs of interest (e.g. with a log(P value) greater than some value such as 10) plus/minus 100 bp.

To pull out genomic intervals, we can use bedtools:
```
module load StdEnv/2020 bedtools/2.29.2
bedtools getfasta [OPTIONS] -fi <input FASTA> -bed <BED/GFF/VCF>
```

# Checking sites with high associatioin with sex
```
cat east_subset_XBwild_unfiltered_allchrs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz.tab | grep 'Chr7S' | egrep '43942|61034|121963|121964|122012|122027|122040|259923|259944|260015|350392|352499|772235|772254|772296|860669|1242454|1242461|1511181|1511264|1511282|1511284|1511316|1531609|1876765|2249073|2249170|2249203|2362031|3263437|3411990|3411991|3919362|3919380|3934882|4627882|4627922|4627936|4713422|4994828|5099318|5099328|5099358|5099415|5887647|6283224|6283268|6283295|6283310|6356456|6356457|6356489|6540061|6593950|7013750|7013758|7234658|7598296|7598300|7709985|7710051|7973851|8022636|8041622|8455125|8552867|8976187|9376010|9448741|9617824|9617833|9635943|9635970|9654309|9654342|10166191|11581482|12242121|12566054|12566110|12590276|12595173|12621554|12766906|13017940|13203777|13203856|14162125|14162166|14374307|14374367|14941293|14981714|14981740|14981816|14981840|15248427|15381479|15617449|15854017|15854025|15854051|15854081|16064130|16064200|16198031|16198038|16198054|16205760|16259300|16355975|16358491|16358587|16684042|16684043|16732802|16732803|17122040|17207897|17367726|17367735|17746179|18122943|18122957|18122959|18852531|20002582|20338155|20366652|20685488|20749096|20756187|20762604|20762642|20762673|20762690|20762711|20762744|20866537|21161007|21161008|22506574|22594345|22895991'
```

# LD
Estimate LD between pairwise comparisons of SNPs:
```
module load nixpkgs/16.09 plink/1.9b_5.2-x86_64
plink --file XB_west_and_Njoro_only_mpileup_raw_wildBorealis_AustinGenome_Chrs_only.vcf.gz_myplink --noweb --r2 dprime --allow-no-sex --allow-extra-chr
```

Get rid of nan lines:
```
cat Njoro_only_mpileup_raw_wildBorealis_AustinGenome_Chrs_only.vcf.gz_myplink.ld | grep -v 'nan' > Njoro_only_mpileup_raw_wildBorealis_AustinGenome_Chrs_only.vcf.gz_myplink.ld_no_nan
```

Plot
```R
dat<-read.table("./XB_west_and_Njoro_only_mpileup_raw_wildBorealis_AustinGenome_Chrs_only.vcf.gz_myplink.ld",header=TRUE)
dat<-read.table("./east_only_mpileup_raw_wildBorealis_AustinGenome_Chrs_only.vcf.gz_myplink_no_nan",header=TRUE)
pdf("./East_only_XLGenome.pdf",w=8, h=12.0, version="1.4", bg="transparent")
p<-ggplot(dat, aes(x=BP_B-BP_A, y=R2)) + 
    # add points
    geom_point(size=2, alpha = 0.7 ) +
    # color the stuff the way I want
    facet_wrap(~CHR_A, ncol = 2) +
    # get rid of gray background
    theme_bw()
p
dev.off()
```
