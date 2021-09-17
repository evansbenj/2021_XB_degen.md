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
