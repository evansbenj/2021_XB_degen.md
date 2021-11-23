# Analysis of genomic windows

Some insights into possible explanations for the mysterious sex chromosomes of XB might be gleaned from analysis of diversity in windows

For this analysis, I used the vcf file that was filtered by GATK but not the one that was then filtered by vcftools to get rid of high coveraeg sites based on a rolling average.  The latter file produced an error with Beagle phasing.

```
#!/bin/sh
#SBATCH --job-name=beagle
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=3:00:00
#SBATCH --mem=128gb
#SBATCH --output=beagle.%J.out
#SBATCH --error=beagle.%J.err
#SBATCH --account=def-ben

# sbatch Beagle.sh chr

module load java

java -Xmx12g -jar /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_dept
h_3sigmas/final_data_including_sites_with_lots_of_missing_data/twisst/beagle.18May20.d20.jar gt=${1} out=${1}_phased.vcf.gz impu
te=true
```

I generated this file using Beagle:
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/WGS_vcfs_by_chr/combined_andfiltereds_gvcfs/allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz_phased.vcf.gz.vcf.gz
```
