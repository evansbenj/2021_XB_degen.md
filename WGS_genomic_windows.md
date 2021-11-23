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

now make a geno file
```
#!/bin/sh                                                                                                      
#SBATCH --job-name=makegeno                                                                                    
#SBATCH --nodes=1                                                                                              
#SBATCH --ntasks-per-node=1                                                                                    
#SBATCH --time=3:00:00                                                                                         
#SBATCH --mem=2gb                                                                                              
#SBATCH --output=makegeno.%J.out                                                                               
#SBATCH --error=makegeno.%J.err                                                                                
#SBATCH --account=def-ben                                                                                      

# sbatch 2020_make_geno_from_vcf.sh path_and_name_of_vcf.gz_file                                               

module load StdEnv/2020
module load scipy-stack/2020b
module load python/3.8.2

python /home/ben/projects/rrg-ben/ben/2017_SEAsian_macaques/SEAsian_macaques_bam/with_papio/2020_Nov_filtered_by_depth_3sigmas/final_data_including_sites_with_lots_of_missing_data/genomics_general/VCF_processing/parseVCF.py -i ${1} | gzip > ${1}.geno.gz
```
