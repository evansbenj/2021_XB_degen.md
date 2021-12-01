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
Calculate pairwise distances
```
./2021_general_genomics_popgen_6pops.sh ../WGS_vcfs_by_chr/combined_andfiltereds_gvcfs/allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz_phased.vcf.gz.vcf.gz.geno.gz East_female East_male West_female West_male Lab_female Lab_male
```

I also did this for chr7S.

I'm plotting these results like this:
```R
library (ggplot2)
library(reshape2) # this facilitates the overlay plot
library(scales) # has function to wrap legends
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2021_Xborealis_sexchr_degen/WGS_genomic_windows")

# read in the data and name the df based on the file name
dat <- read.table("East_female_East_male_West_female_West_male_Lab_female_Lab_male_windowstats_100000bp.csv", header = T, sep=",")

new_dat <- dat[,c(1,2,6:11)]
new_dat$East <- new_dat$pi_East_female/new_dat$pi_East_male
new_dat$Lab <- new_dat$pi_Lab_female/new_dat$pi_Lab_male
new_dat$West <- new_dat$pi_West_female/new_dat$pi_West_male
new_dat$SL_or_not <- "light blue" # not_SL
new_dat$SL_or_not[(new_dat$scaffold == "Chr8L") & (new_dat$start < 54000000)] <- "red" # this is SL
#View(new_dat)
# subset the data to include only the position, ratios, and SL column
subset_dat <- new_dat[,c(1,2,9:12)]
# need to melt the data
melted <- melt(subset_dat, id=c("scaffold","start","SL_or_not"))

# ok now determine the order that we want the variables to be printed
melted$SL_or_not <- factor(melted$SL_or_not,
                       levels = c("red","light blue"),ordered = TRUE)
melted$variable <- factor(melted$variable,
                           levels = c("East",
                                      "Lab",
                                      "West"),ordered = TRUE)
# Now make a box and whisker plot
# Notched Boxplot of Tooth Growth Against 2 Crossed Factors
# boxes colored for ease of interpretation
ggplot(melted, aes(x=variable, y=value, fill=SL_or_not)) + 
  geom_boxplot() +
  theme_bw() +
  scale_fill_manual(values = c("red", "light blue"),
                      name="<54Mb?",
                      breaks=c("red", "light blue"),
                      labels=c("Yes", "No")) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="F:M polymorphism ratio")
```
