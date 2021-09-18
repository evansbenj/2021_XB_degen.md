# Plink

Associations between SNPs and a phenotype (such as phenotypic sex) are easily calculated using plink.

```
module load plink/1.9b_5.2-x86_64
module load StdEnv/2020
module load r/4.0.2
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
``
module load StdEnv/2020  gcc/9.3.0 bcftools/1.11

bcftools view mpileup_raw_wildBorealis_AuXXXGenome.vcf.gz --regions Chr1L,Chr1S,Chr2L,Chr2S,Chr3L,Chr3S,Chr4L,Chr4S,Chr5L,Chr5S,Chr6L,Chr6S,Chr7L,Chr7S,Chr8L,Chr8S,Chr9_10L,Chr9_10S -O z -o mpileup_raw_wildBorealis_AuXXXGenome_Chrs_only.vcf.gz
``

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
bcftools view -s Fem_Cheplaskei_BJE4459_Xb.fq.gz,Fem_Cheplaskei_BJE4461_Xb.fq.gz,Fem_Cheplaskei_BJE4470_Xb.fq.gz,Fem_Chesuwe_BJE4479_Xb.fq.gz,Fem_Chesuwe_BJE4481_Xb.fq.gz,Fem_Eldoret_BJE4471_Xb.fq.gz,Fem_Eldoret_BJE4472_Xb.fq.gz,Fem_Eldoret_BJE4474_Xb.fq.gz,Fem_Eldoret_BJE4475_Xb.fq.gz,Fem_Eldoret_BJE4476_Xb.fq.gz,Fem_Kiminini_BJE4429_Xb.fq.gz,Fem_Kiminini_BJE4433_Xb.fq.gz,Fem_Lukhome_BJE4441_Xb.fq.gz,Fem_Lukhome_BJE4444_Xb.fq.gz,Fem_Lukhome_BJE4445_Xb.fq.gz,Fem_Lukhome_BJE4446_Xb.fq.gz,Mal_Cheplaskei_BJE4460_Xb.fq.gz,Mal_Cheplaskei_BJE4462_Xb.fq.gz,Mal_Cheplaskei_BJE4465_Xb.fq.gz,Mal_Cheplaskei_BJE4469_Xb.fq.gz,Mal_Chesuwe_BJE4477_Xb.fq.gz,Mal_Chesuwe_BJE4478_Xb.fq.gz,Mal_Chesuwe_BJE4480_Xb.fq.gz,Mal_Eldoret_BJE4473_Xb.fq.gz,Mal_Kiminini_BJE4430_Xb.fq.gz,Mal_Kiminini_BJE4431_Xb.fq.gz,Mal_Kiminini_BJE4432_Xb.fq.gz,Mal_Kisumu_BJE4391_Xb.fq.gz,Mal_Lukhome_BJE4442_Xb.fq.gz,Mal_Lukhome_BJE4443_Xb.fq.gz,Mal_Lukhome_BJE4447_Xb.fq.gz mpileup_raw_wildBorealis_AuXXXnGenome_Chrs_only.vcf.gz -O z > west_only_mpileup_raw_wildBorealis_AuXXXGenome_Chrs_only.vcf.gz

bcftools view -s Fem_NMK_BJE4562_Xb.fq.gz,Fem_NMK_BJE4564_Xb.fq.gz,Fem_NMK_BJE4567_Xb.fq.gz,Fem_NMK_BJE4568_Xb.fq.gz,Fem_Nakuru_BJE4363_Xb.fq.gz,Fem_Nakuru_BJE4364_Xb.fq.gz,Fem_Nakuru_BJE4367_Xb.fq.gz,Fem_Nakuru_BJE4368_Xb.fq.gz,Mal_NMK_BJE4563_Xb.fq.gz,Mal_Nakuru_BJE4365_Xb.fq.gz,Mal_Nakuru_BJE4366_Xb.fq.gz,Mal_Nakuru_BJE4374_Xb.fq.gz,Mal_Nakuru_BJE4377_Xb.fq.gz mpileup_raw_wildBorealis_AuXXXGenome_Chrs_only.vcf.gz -O z > central_only_mpileup_raw_wildBorealis_AuXXXGenome_Chrs_only.vcf.gz

bcftools view -s Fem_Wundanyi_BJE4515_Xb.fq.gz,Fem_Wundanyi_BJE4516_Xb.fq.gz,Fem_Wundanyi_BJE4534_Xb.fq.gz,Fem_Wundanyi_BJE4535_Xb.fq.gz,Fem_Wundanyi_BJE4541_Xb.fq.gz,Mal_Wundanyi_BJE4536_Xb.fq.gz,Mal_Wundanyi_BJE4537_Xb.fq.gz,Mal_Wundanyi_BJE4538_Xb.fq.gz,Mal_Wundanyi_BJE4539_Xb.fq.gz,Mal_Wundanyi_BJE4540_Xb.fq.gz mpileup_raw_wildBorealis_AuXXXGenome_Chrs_only.vcf.gz -O z > east_wund_only_mpileup_raw_wildBorealis_AuXXXGenome_Chrs_only.vcf.gz

bcftools view -s Fem_Nakuru_BJE4363_Xb.fq.gz,Fem_Nakuru_BJE4364_Xb.fq.gz,Fem_Nakuru_BJE4367_Xb.fq.gz,Fem_Nakuru_BJE4368_Xb.fq.gz,Mal_Nakuru_BJE4365_Xb.fq.gz,Mal_Nakuru_BJE4366_Xb.fq.gz,Mal_Nakuru_BJE4374_Xb.fq.gz,Mal_Nakuru_BJE4377_Xb.fq.gz mpileup_raw_wildBorealis_AuXXXGenome_Chrs_only.vcf.gz -O z > Njoror_only_mpileup_raw_wildBorealis_AuXXXGenome_Chrs_only.vcf.gz

bcftools view -s Fem_Cheplaskei_BJE4459_Xb.fq.gz,Fem_Cheplaskei_BJE4461_Xb.fq.gz,Fem_Cheplaskei_BJE4470_Xb.fq.gz,Fem_Chesuwe_BJE4479_Xb.fq.gz,Fem_Chesuwe_BJE4481_Xb.fq.gz,Fem_Eldoret_BJE4471_Xb.fq.gz,Fem_Eldoret_BJE4472_Xb.fq.gz,Fem_Eldoret_BJE4474_Xb.fq.gz,Fem_Eldoret_BJE4475_Xb.fq.gz,Fem_Eldoret_BJE4476_Xb.fq.gz,Fem_Kiminini_BJE4429_Xb.fq.gz,Fem_Kiminini_BJE4433_Xb.fq.gz,Fem_Lukhome_BJE4441_Xb.fq.gz,Fem_Lukhome_BJE4444_Xb.fq.gz,Fem_Lukhome_BJE4445_Xb.fq.gz,Fem_Lukhome_BJE4446_Xb.fq.gz,Fem_Nakuru_BJE4363_Xb.fq.gz,Fem_Nakuru_BJE4364_Xb.fq.gz,Fem_Nakuru_BJE4367_Xb.fq.gz,Fem_Nakuru_BJE4368_Xb.fq.gz,Mal_Cheplaskei_BJE4460_Xb.fq.gz,Mal_Cheplaskei_BJE4462_Xb.fq.gz,Mal_Cheplaskei_BJE4465_Xb.fq.gz,Mal_Cheplaskei_BJE4469_Xb.fq.gz,Mal_Chesuwe_BJE4477_Xb.fq.gz,Mal_Chesuwe_BJE4478_Xb.fq.gz,Mal_Chesuwe_BJE4480_Xb.fq.gz,Mal_Eldoret_BJE4473_Xb.fq.gz,Mal_Kiminini_BJE4430_Xb.fq.gz,Mal_Kiminini_BJE4431_Xb.fq.gz,Mal_Kiminini_BJE4432_Xb.fq.gz,Mal_Kisumu_BJE4391_Xb.fq.gz,Mal_Lukhome_BJE4442_Xb.fq.gz,Mal_Lukhome_BJE4443_Xb.fq.gz,Mal_Lukhome_BJE4447_Xb.fq.gz,Mal_Nakuru_BJE4365_Xb.fq.gz,Mal_Nakuru_BJE4366_Xb.fq.gz,Mal_Nakuru_BJE4374_Xb.fq.gz,Mal_Nakuru_BJE4377_Xb.fq.gz mpileup_raw_wildBorealis_AuXXXnGenome_Chrs_only.vcf.gz -O z > XB_west_and_Njoro_only_mpileup_raw_wildBorealis_AuXXXGenome_Chrs_only.vcf.gz

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
