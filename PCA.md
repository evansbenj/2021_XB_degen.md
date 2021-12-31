Working directory:
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/XB_RADseq_mapped_to_AustinXB_BJE_2021/combined_and_genotyped_vcfs_trim
```
subset chr7S SL and nonSL:
```
vcftools --gzvcf XBwild_unfiltered_allchrs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz --chr Chr7S --from-bp 1 --to-bp 23000000 --out XBwild_unfiltered_allchrs.vcf.gz_filtered.vcf.gz_filtered_removed_Chr7S_SL --recode
```
```
vcftools --gzvcf XBwild_unfiltered_allchrs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz --chr Chr7S --from-bp 23000001 --to-bp 105895007 --out XBwild_unfiltered_allchrs.vcf.gz_filtered.vcf.gz_filtered_removed_Chr7S_notSL --recode
```
Do PCA analysis
```R

```
