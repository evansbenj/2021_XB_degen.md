# Making a nj tree for SL and non-SL portion of Chr8L from WGS data
Directory:
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/WGS_vcfs_by_chr/combined_andfiltereds_gvcfs
```

First make vcf files with SL and non-SL portions of Chr8L from WGS data:
```
bcftools view allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz --regions Chr8L:1-54100000 > allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_SLonly.vcf
bcftools view allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz --regions Chr8L:54100001-123836260 > allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_pseudoonly.vcf
```
Now make a nexus file (-n) from each one, and require that everyone have a genotype:
```
python /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/laevis/bin/vcf2phylip/vcf2phylip.py -i allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_SLonly.vcf -n --min-samples-locus 6 --output-prefix Chr8L_SLonly_nexus
python /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_allo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/laevis/bin/vcf2phylip/vcf2phylip.py -i allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_pseudoonly.vcf -n --min-samples-locus 6 --output-prefix Chr8L_nonSLonly_nexus
```

Resulting files:
```
Chr8L_SLonly_nexus.min6.nex
Chr8L_nonSLonly_nexus.min6.nex
```
