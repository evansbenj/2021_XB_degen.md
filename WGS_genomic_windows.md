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

# More plots Dec 23, 2021
```R
library (ggplot2)
library(reshape2) # this facilitates the overlay plot
library(scales) # has function to wrap legends
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2021_Xborealis_sexchr_degen/WGS_genomic_windows")

# read in the data and name the df based on the file name
dat <- read.table("East_female_East_male_West_female_West_male_Lab_female_Lab_male_windowstats_100000bp_Chr8L.csv", header = T, sep=",")

new_dat <- dat[,c(1,3,6:11)]
new_dat$East <- new_dat$pi_East_female/new_dat$pi_East_male
new_dat$Lab <- new_dat$pi_Lab_female/new_dat$pi_Lab_male
new_dat$West <- new_dat$pi_West_female/new_dat$pi_West_male
new_dat$SL_or_not <- "light blue" # not_SL
new_dat$SL_or_not[(new_dat$scaffold == "Chr8L") & (new_dat$end < 54100000)] <- "red" # this is SL
#View(new_dat)
# subset the data to include only the position, ratios, and SL column
subset_dat <- new_dat[,c(1,2,9:12)]
# need to melt the data
melted <- melt(subset_dat, id=c("scaffold","end","SL_or_not"))

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
png("F_to_M_PolymorphismRatios_Chr8L.png", width = 110, height = 110, units='mm', res = 100) 
  ggplot(melted, aes(x=variable, y=log(value), fill=SL_or_not)) + 
  geom_hline(yintercept=0, linetype="solid",color = "black", size=0.5) + 
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("red", "light blue"),
                      name="<54.1 Mb?",
                      breaks=c("red", "light blue"),
                      labels=c("Yes", "No")) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="log(F:M polymorphism ratio)") 
dev.off()
  

# read in the data and name the df based on the file name
dat <- read.table("East_female_East_male_West_female_West_male_Lab_female_Lab_male_windowstats_100000bp_Chr7S.csv", header = T, sep=",")

new_dat <- dat[,c(1,3,6:11)]
new_dat$East <- new_dat$pi_East_female/new_dat$pi_East_male
new_dat$Lab <- new_dat$pi_Lab_female/new_dat$pi_Lab_male
new_dat$West <- new_dat$pi_West_female/new_dat$pi_West_male
new_dat$SL_or_not <- "light blue" # not_SL
new_dat$SL_or_not[(new_dat$scaffold == "Chr7S") & (new_dat$end < 23000000)] <- "red" # this is SL
#View(new_dat)
# subset the data to include only the position, ratios, and SL column
subset_dat <- new_dat[,c(1,2,9:12)]
# need to melt the data
melted <- melt(subset_dat, id=c("scaffold","end","SL_or_not"))

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
png("F_to_M_PolymorphismRatios_Chr7S.png", width = 110, height = 110, units='mm', res = 100)  
ggplot(melted, aes(x=variable, y=log(value), fill=SL_or_not)) + 
  geom_hline(yintercept=0, linetype="solid",color = "black", size=0.5) + 
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("red", "light blue"),
                    name="<23Mb?",
                    breaks=c("red", "light blue"),
                    labels=c("Yes", "No")) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="log(F:M polymorphism ratio)")
dev.off()


# Now plot the Male_Lab:Male_east polymorphism ratio
rm(list = ls())
dat <- read.table("East_female_East_male_West_female_West_male_Lab_female_Lab_male_windowstats_100000bp_Chr8L.csv", header = T, sep=",")
new_dat <- dat[,c(1,3,6:11)]
new_dat$M_Lab_to_M_east <- new_dat$pi_Lab_male/new_dat$pi_East_male

new_dat$SL_or_not <- "light blue" # not_SL
new_dat$SL_or_not[(new_dat$scaffold == "Chr8L") & (new_dat$end < 54000000)] <- "red" # this is SL
#View(new_dat)
# subset the data to include only the position, ratios, and SL column
subset_dat <- new_dat[,c(1,2,9:10)]
# need to melt the data
melted <- melt(subset_dat, id=c("scaffold","end","SL_or_not"))

# ok now determine the order that we want the variables to be printed
melted$SL_or_not <- factor(melted$SL_or_not,
                           levels = c("red","light blue"),ordered = TRUE)

# boxes colored for ease of interpretation
png(filename = "MLab_to_Meast_PolymorphismRatios_Chr8L.png",w=300, h=300,units = "px", bg="transparent")
ggplot(melted, aes(x=variable, y=value, fill=SL_or_not)) + 
  geom_boxplot() +
  theme_classic() +
  scale_fill_manual(values = c("red", "light blue"),
                    name="<54Mb?",
                    breaks=c("red", "light blue"),
                    labels=c("Yes", "No")) +
  scale_x_discrete(name="") +
  scale_y_continuous(name="F:M polymorphism ratio")
dev.off()


# Now make boxplots of pi by individual in SL and nonSL regions for each chr
# Chr8L
# clear the environment
rm(list = ls())

# read in the data and name the df based on the file name
dat <- read.table("East_female_East_male_West_female_West_male_Lab_female_Lab_male_windowstats_100000bp_Chr8L.csv", header = T, sep=",")
# subset each dataset
west_female <- dat[,c(1,3,8)]
west_male <- dat[,c(1,3,9)]
east_female <- dat[,c(1,3,6)]
east_male <- dat[,c(1,3,7)]
lab_male <- dat[,c(1,3,11)]
lab_female <- dat[,c(1,3,10)]
# add a sample variable
west_female$Individual <- "west female"
west_male$Individual <- "west male"
east_female$Individual <- "east female"
east_male$Individual <- "east male"
lab_female$Individual <- "lab female"
lab_male$Individual <- "lab male"
# fix the column names
colnames(lab_male)[3] <- "Pi"
colnames(lab_female)[3] <- "Pi"
colnames(west_male)[3] <- "Pi"
colnames(west_female)[3] <- "Pi"
colnames(east_male)[3] <- "Pi"
colnames(east_female)[3] <- "Pi"
# combine them
newdat <- rbind(east_female,east_male,lab_female,lab_male,west_female,west_male)
# make a column for SL_or_not
newdat$SL_or_not <- "light blue" # not_SL
newdat$SL_or_not[(newdat$scaffold == "Chr8L") & (newdat$end < 54100000)] <- "red" # this is SL
# get some numbers (median depth of each section for each sample)


#east_female_SL
median(newdat$Pi[(newdat$Individual == "East female") & 
                      (newdat$scaffold == "Chr8L") 
                    & (newdat$end < 54100000)], na.rm=T)
# [1] 0.3648
#east_female_notSL
median(newdat$Pi[(newdat$Individual == "East female") & 
                   (newdat$scaffold == "Chr8L") 
                 & (newdat$end >= 54100000)], na.rm=T)
# [1] 0.1186
#east_male_SL
median(newdat$Pi[(newdat$Individual == "East male") & 
                   (newdat$scaffold == "Chr8L") 
                 & (newdat$end < 54100000)], na.rm=T)
# [1] 0.11745
#east_male_notSL
median(newdat$Pi[(newdat$Individual == "East male") & 
                   (newdat$scaffold == "Chr8L") 
                 & (newdat$end >= 54100000)], na.rm=T)
# [1] 0.1192
#west_female_SL
median(newdat$Pi[(newdat$Individual == "West female") & 
                   (newdat$scaffold == "Chr8L") 
                 & (newdat$end < 54100000)], na.rm=T)
# [1] 0.24485
#west_female_notSL
median(newdat$Pi[(newdat$Individual == "West female") & 
                   (newdat$scaffold == "Chr8L") 
                 & (newdat$end >= 54100000)], na.rm=T)
# [1] 0.2422
#west_female_SL
median(newdat$Pi[(newdat$Individual == "West male") & 
                   (newdat$scaffold == "Chr8L") 
                 & (newdat$end < 54100000)], na.rm=T)
# [1] 0.23165
#west_female_notSL
median(newdat$Pi[(newdat$Individual == "West male") & 
                   (newdat$scaffold == "Chr8L") 
                 & (newdat$end >= 54100000)], na.rm=T)
# [1] 0.2482


# make the SL first
newdat$SL_or_not <- factor(newdat$SL_or_not,
                           levels = c("red","light blue"),ordered = TRUE)

png("SL_notSL_pi_by_individual_Chr8L.png", width = 200, height = 110, units='mm', res = 100) 
ggplot(newdat, aes(x=Individual, y=Pi, fill=SL_or_not)) + 
  geom_boxplot(outlier.size = 0.5) +
  #coord_cartesian(ylim = c(0, 100)) +
  theme_classic() +
  scale_fill_manual(values = c("red", "light blue"),
                    name="<54.1 Mb?",
                    breaks=c("red", "light blue"),
                    labels=c("Yes", "No")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  scale_y_continuous(name="Pairwise nucleotide diversity (Chr8L)")
dev.off()



# Now make boxplots of pi by individual in SL and nonSL regions for each chr
# Chr7S
# clear the environment
rm(list = ls())

# read in the data and name the df based on the file name
dat <- read.table("East_female_East_male_West_female_West_male_Lab_female_Lab_male_windowstats_100000bp_Chr7S.csv", header = T, sep=",")
# subset each dataset
west_female <- dat[,c(1,3,8)]
west_male <- dat[,c(1,3,9)]
east_female <- dat[,c(1,3,6)]
east_male <- dat[,c(1,3,7)]
lab_male <- dat[,c(1,3,11)]
lab_female <- dat[,c(1,3,10)]
# add a sample variable
west_female$Individual <- "West female"
west_male$Individual <- "West male"
east_female$Individual <- "East female"
east_male$Individual <- "East male"
lab_female$Individual <- "Lab female"
lab_male$Individual <- "Lab male"
# fix the column names
colnames(lab_male)[3] <- "Pi"
colnames(lab_female)[3] <- "Pi"
colnames(west_male)[3] <- "Pi"
colnames(west_female)[3] <- "Pi"
colnames(east_male)[3] <- "Pi"
colnames(east_female)[3] <- "Pi"
# combine them
newdat <- rbind(east_female,east_male,lab_female,lab_male,west_female,west_male)
# make a column for SL_or_not
newdat$SL_or_not <- "light blue" # not_SL
newdat$SL_or_not[(newdat$scaffold == "Chr7S") & (newdat$end < 23000000)] <- "red" # this is SL
# get some numbers (median depth of each section for each sample)


#east_female_SL
median(newdat$Pi[(newdat$Individual == "East female") & 
                   (newdat$scaffold == "Chr7S") 
                 & (newdat$end < 23000000)], na.rm=T)
# [1] 0.1165
#east_female_notSL
median(newdat$Pi[(newdat$Individual == "East female") & 
                   (newdat$scaffold == "Chr7S") 
                 & (newdat$end >= 23000000)], na.rm=T)
# [1] 0.1235
#east_male_SL
median(newdat$Pi[(newdat$Individual == "East male") & 
                   (newdat$scaffold == "Chr7S") 
                 & (newdat$end < 23000000)], na.rm=T)
# [1] 0.3536
#east_male_notSL
median(newdat$Pi[(newdat$Individual == "East male") & 
                   (newdat$scaffold == "Chr7S") 
                 & (newdat$end >= 23000000)], na.rm=T)
# [1] 0.12325
#west_female_SL
median(newdat$Pi[(newdat$Individual == "West female") & 
                   (newdat$scaffold == "Chr7S") 
                 & (newdat$end < 23000000)], na.rm=T)
# [1] 0.2301
#west_female_notSL
median(newdat$Pi[(newdat$Individual == "West female") & 
                   (newdat$scaffold == "Chr7S") 
                 & (newdat$end >= 23000000)], na.rm=T)
# [1] 0.25595
#west_female_SL
median(newdat$Pi[(newdat$Individual == "West male") & 
                   (newdat$scaffold == "Chr7S") 
                 & (newdat$end < 23000000)], na.rm=T)
# [1] 0.2414
#west_female_notSL
median(newdat$Pi[(newdat$Individual == "West male") & 
                   (newdat$scaffold == "Chr7S") 
                 & (newdat$end >= 23000000)], na.rm=T)
# [1] 0.2563


# make the SL first
newdat$SL_or_not <- factor(newdat$SL_or_not,
                           levels = c("red","light blue"),ordered = TRUE)

png("SL_notSL_pi_by_individual_Chr7S.png", width = 200, height = 110, units='mm', res = 100) 
ggplot(newdat, aes(x=Individual, y=Pi, fill=SL_or_not)) + 
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("red", "light blue"),
                    name="<23.0 Mb?",
                    breaks=c("red", "light blue"),
                    labels=c("Yes", "No")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  scale_y_continuous(name="Pairwise nucleotide diversity (Chr7S)")
dev.off()


# now make barplot from RRGS data
# Chr7S
# clear the environment
rm(list = ls())

# read in the data and name the df based on the file name
dat <- read.table("Fem_BJE4515_Fem_BJE4516_Fem_BJE4534_Fem_BJE4535_Fem_BJE4541_Mal_BJE4536_Mal_BJE4537_Mal_BJE4538_Mal_BJE4539_Mal_BJE4540_windowstats_2000000bp_Chr7S.csv", header = T, sep=",")
# subset each dataset
Fem_BJE4515 <- dat[,c(1,3,6)]
Fem_BJE4516 <- dat[,c(1,3,7)]
Fem_BJE4534 <- dat[,c(1,3,8)]
Fem_BJE4535 <- dat[,c(1,3,9)]
Fem_BJE4541 <- dat[,c(1,3,10)]
Mal_BJE4536 <- dat[,c(1,3,11)]
Mal_BJE4537 <- dat[,c(1,3,12)]
Mal_BJE4538 <- dat[,c(1,3,13)]
Mal_BJE4539 <- dat[,c(1,3,14)]
Mal_BJE4540 <- dat[,c(1,3,15)]
# add a sample variable
Fem_BJE4515$Individual <- "BJE4515"
Fem_BJE4516$Individual <- "BJE4516"
Fem_BJE4534$Individual <- "BJE4534"
Fem_BJE4535$Individual <- "BJE4535"
Fem_BJE4541$Individual <- "BJE4541"
Mal_BJE4536$Individual <- "BJE4536"
Mal_BJE4537$Individual <- "BJE4537"
Mal_BJE4538$Individual <- "BJE4538"
Mal_BJE4539$Individual <- "BJE4539"
Mal_BJE4540$Individual <- "BJE4540"
# add a sex variable
Fem_BJE4515$Sex <- "Females"
Fem_BJE4516$Sex <- "Females"
Fem_BJE4534$Sex <- "Females"
Fem_BJE4535$Sex <- "Females"
Fem_BJE4541$Sex <- "Females"
Mal_BJE4536$Sex <- "Males"
Mal_BJE4537$Sex <- "Males"
Mal_BJE4538$Sex <- "Males"
Mal_BJE4539$Sex <- "Males"
Mal_BJE4540$Sex <- "Males"
# fix the column names
colnames(Fem_BJE4515)[3] <- "Pi"
colnames(Fem_BJE4516)[3] <- "Pi"
colnames(Fem_BJE4534)[3] <- "Pi"
colnames(Fem_BJE4535)[3] <- "Pi"
colnames(Fem_BJE4541)[3] <- "Pi"
colnames(Mal_BJE4536)[3] <- "Pi"
colnames(Mal_BJE4537)[3] <- "Pi"
colnames(Mal_BJE4538)[3] <- "Pi"
colnames(Mal_BJE4539)[3] <- "Pi"
colnames(Mal_BJE4540)[3] <- "Pi"

# combine them
newdat <- rbind(Fem_BJE4515,Fem_BJE4516,Fem_BJE4534,Fem_BJE4535,Fem_BJE4541,Mal_BJE4536,Mal_BJE4537,
                Mal_BJE4538,Mal_BJE4539,Mal_BJE4540)
# make a column for SL_or_not
newdat$SL_or_not <- "light blue" # not_SL
newdat$SL_or_not[(newdat$scaffold == "Chr7S") & (newdat$end < 23000000)] <- "red" # this is SL


# make the SL first
newdat$SL_or_not <- factor(newdat$SL_or_not,
                           levels = c("red","light blue"),ordered = TRUE)

png("SL_notSL_pi_by_individual_RRGS_Chr7S.png", width = 200, height = 110, units='mm', res = 100) 
ggplot(newdat, aes(x=Individual, y=Pi, fill=SL_or_not)) + 
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("red", "light blue"),
                    name="<23.0 Mb?",
                    breaks=c("red", "light blue"),
                    labels=c("Yes", "No")) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14)) +
  scale_y_continuous(name="Pairwise nucleotide diversity (Chr7S)") +
  facet_wrap(. ~ Sex, scales = "free_x") +
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(size = 14))+
  theme(legend.position = c(0.95, 0.8))
dev.off()


# now make barplot from RRGS data
# Chr8L
# clear the environment
rm(list = ls())

# read in the data and name the df based on the file name
dat <- read.table("Fem_BJE4515_Fem_BJE4516_Fem_BJE4534_Fem_BJE4535_Fem_BJE4541_Mal_BJE4536_Mal_BJE4537_Mal_BJE4538_Mal_BJE4539_Mal_BJE4540_windowstats_2000000bp_Chr8L.csv", header = T, sep=",")
# subset each dataset
Fem_BJE4515 <- dat[,c(1,3,6)]
Fem_BJE4516 <- dat[,c(1,3,7)]
Fem_BJE4534 <- dat[,c(1,3,8)]
Fem_BJE4535 <- dat[,c(1,3,9)]
Fem_BJE4541 <- dat[,c(1,3,10)]
Mal_BJE4536 <- dat[,c(1,3,11)]
Mal_BJE4537 <- dat[,c(1,3,12)]
Mal_BJE4538 <- dat[,c(1,3,13)]
Mal_BJE4539 <- dat[,c(1,3,14)]
Mal_BJE4540 <- dat[,c(1,3,15)]
# add a sample variable
Fem_BJE4515$Individual <- "BJE4515"
Fem_BJE4516$Individual <- "BJE4516"
Fem_BJE4534$Individual <- "BJE4534"
Fem_BJE4535$Individual <- "BJE4535"
Fem_BJE4541$Individual <- "BJE4541"
Mal_BJE4536$Individual <- "BJE4536"
Mal_BJE4537$Individual <- "BJE4537"
Mal_BJE4538$Individual <- "BJE4538"
Mal_BJE4539$Individual <- "BJE4539"
Mal_BJE4540$Individual <- "BJE4540"
# add a sex variable
Fem_BJE4515$Sex <- "Females"
Fem_BJE4516$Sex <- "Females"
Fem_BJE4534$Sex <- "Females"
Fem_BJE4535$Sex <- "Females"
Fem_BJE4541$Sex <- "Females"
Mal_BJE4536$Sex <- "Males"
Mal_BJE4537$Sex <- "Males"
Mal_BJE4538$Sex <- "Males"
Mal_BJE4539$Sex <- "Males"
Mal_BJE4540$Sex <- "Males"
# fix the column names
colnames(Fem_BJE4515)[3] <- "Pi"
colnames(Fem_BJE4516)[3] <- "Pi"
colnames(Fem_BJE4534)[3] <- "Pi"
colnames(Fem_BJE4535)[3] <- "Pi"
colnames(Fem_BJE4541)[3] <- "Pi"
colnames(Mal_BJE4536)[3] <- "Pi"
colnames(Mal_BJE4537)[3] <- "Pi"
colnames(Mal_BJE4538)[3] <- "Pi"
colnames(Mal_BJE4539)[3] <- "Pi"
colnames(Mal_BJE4540)[3] <- "Pi"

# combine them
newdat <- rbind(Fem_BJE4515,Fem_BJE4516,Fem_BJE4534,Fem_BJE4535,Fem_BJE4541,Mal_BJE4536,Mal_BJE4537,
                Mal_BJE4538,Mal_BJE4539,Mal_BJE4540)
# make a column for SL_or_not
newdat$SL_or_not <- "light blue" # not_SL
newdat$SL_or_not[(newdat$scaffold == "Chr8L") & (newdat$end < 54100000)] <- "red" # this is SL


# make the SL first
newdat$SL_or_not <- factor(newdat$SL_or_not,
                           levels = c("red","light blue"),ordered = TRUE)

png("SL_notSL_pi_by_individual_RRGS_Chr8L.png", width = 200, height = 110, units='mm', res = 100) 
ggplot(newdat, aes(x=Individual, y=Pi, fill=SL_or_not)) + 
  geom_boxplot(outlier.size = 0.5) +
  theme_classic() +
  scale_fill_manual(values = c("red", "light blue"),
                    name="<54.1 Mb?",
                    breaks=c("red", "light blue"),
                    labels=c("Yes", "No")) +
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=14)) +
  scale_y_continuous(name="Pairwise nucleotide diversity (Chr8L)") +
  facet_wrap(. ~ Sex, scales = "free_x") +
  theme(strip.background = element_blank())+
  theme(strip.text.x = element_text(size = 14)) +
  theme(legend.position = c(0.95, 0.8)) #+
  #annotate("segment",x = 0.5,xend = 5.5,y = 0.23,yend = 0.23, size = 0.5) 
 # geom_segment(aes(x =0.5, y = 0.23, xend = 5.5, yend = 0.23), 
#               colour = "black", size = .5,
 #              data = newdat)
dev.off()
```
```R
library (ggplot2)
library(reshape2) # this facilitates the overlay plot
library(scales) # has function to wrap legends
setwd("/Users/Shared/Previously Relocated Items/Security/projects/2021_Xborealis_sexchr_degen/WGS_genomic_windows")

# This is the heterozygous count data from my perl script that counts sex-limited and sex-shared
# polymorphisms in the west genotypes; below is code that plots genomic windows for Chr8L
# read in the data and name the df based on the file name
XB_East_West_windows <- read.table("EastWest_windows_counts_lt54MB.txt", header = F, sep="\t")
colnames(XB_East_West_windows) <- c("Lower","Upper","Individual",
                                    "fem_n_homoz_female_specific","fem_n_homoz_sex_shared",
                                    "fem_n_het","Individual","mal_n_homoz_female_specific",
                                    "mal_n_homoz_sex_shared","mal_n_het","Conserved","Control")
female_data <-XB_East_West_windows[,c(1:6,12)]
male_data <-XB_East_West_windows[,c(1:2,7:10,14)]
colnames(female_data) <- c("Lower","Upper","Individual",
                                    "n_homoz_female_specific","n_homoz_sex_shared",
                                    "n_het","Conserved_Control")
colnames(male_data) <- c("Lower","Upper","Individual",
                           "n_homoz_female_specific","n_homoz_sex_shared",
                           "n_het","Conserved_Control") 

EastWest_windows_plot_data <- rbind(female_data,male_data)
# The interesting ones are the ones that are heterozygous in all three females from 
# all three locations. This would be consistent with no change in the sex-determining
# gene but a big change in recombination suppression
# So I calculated the three genotypes in the west female and male for the
# female-specific and sex-shared SNPs that are in both females (East and Lab, for the former)
# and (for the sex shared) in both females and homoz in both males (East and Lab).

# the expectation is that if the sex determining locus is the same in all populations
# that there should be a region where all females and no males are heterozygous.

# measuring this in 100,000 bp windows, which is pretty small

#png(filename = "EastWest_windows_female_specific_hets.png",w=1200, h=400,units = "px", bg="transparent")
pdf("./EastWest_windows_female_specific_hets.pdf",w=8, h=4.0, version="1.4", bg="transparent")
p<-ggplot(EastWest_windows_plot_data, aes(x=Lower/1000000, y=Conserved_Control, col = Individual)) + 
  #geom_point(aes(color = factor(variable))) +
  scale_color_manual(labels = c("West female",
                                "West male"),  
  values=c("red","blue")) +
  geom_point(size=0.25) + #, alpha = 0.2) +
  geom_line()+
  #geom_line(aes(colour = "red"), linetype = 1) +
  scale_y_continuous(name="Count per 100,000 bp window") + #, limits=c(-100,100)) +
  # make that legend wrap
  guides(linetype=guide_legend(nrow=2)) +
  # log transform y-axis
  scale_x_continuous(name="Position on Chr8L (Mb)", 
                     breaks=seq(0,57,10), limits = c(1,54)) +
  # add a box to highlight the cool section
  #geom_rect(data=EastWest_windows_plot_data, mapping=aes(xmin=c(26.7), xmax=c(30), 
  #          ymin=c(0), ymax=c(60)), fill=NA,
  #          color="red", alpha=0.5) +
  geom_segment(aes(x = 27.207878, y = 56, xend = 27.207878, yend = 60), colour = "black", 
               data = EastWest_windows_plot_data) + 
  geom_text(aes(x=26, y= 63.4, label = "mtmr8", angle=-45), colour = "black", vjust = 1, size = 3) +
  geom_segment(aes(x = 28.626595, y = 56, xend = 28.626595, yend = 60), colour = "black",
               data = EastWest_windows_plot_data) + 
  geom_text(aes(x=27.6, y= 63.4, label = "abcb7", angle=-45), colour = "black", vjust = 1, size = 3) +
  #geom_vline(data = data.frame(xint=28.626595), aes(xintercept = xint), linetype = "solid", lwd=1) + #abcb7
  geom_segment(aes(x =29.266952, y = 56, xend = 29.266952, yend = 60), colour = "black", 
               data = EastWest_windows_plot_data) + 
  geom_segment(aes(x =29.278133, y = 56, xend = 29.278133, yend = 60), colour = "black", 
               data = EastWest_windows_plot_data) + 
  geom_text(aes(x=29, y= 64, label = "reg1/2", angle=-45), colour = "black", vjust = 1, size = 3) +
  geom_segment(aes(x =37.850413, y = 56, xend = 37.850413, yend = 60), colour = "black", 
               data = EastWest_windows_plot_data) + 
  geom_text(aes(x=37.6, y= 63.4, label = "sox3", angle=-45), colour = "black", vjust = 1, size = 3) +
  geom_segment(aes(x =12.084490, y = 56, xend = 12.084490, yend = 60), colour = "black", 
               data = EastWest_windows_plot_data) + 
  geom_text(aes(x=11.5, y= 63.4, label = "nr5a1", angle=-45), colour = "black", vjust = 1, size = 3) +
  # get rid of gray background
  theme_classic() +
  theme(legend.position = c(0.9, 0.8))
  #theme(legend.position = "none")
  #facet_wrap( ~ pop, nrow=4) +theme(strip.text.x = element_blank())
p
dev.off()




# Get max window
EastWest_windows_plot_data[EastWest_windows_plot_data$Conserved_Control==max(EastWest_windows_plot_data$Conserved_Control),]




# read in the data and name the df based on the file name
XB_East_West <- read.table("East_female_East_male_West_female_West_male_Lab_female_Lab_male_windowstats_100000bp_Chr8L.csv", header = T, sep=",")

##### PI #####

# make a new df of pi for overlay plots
pi_East_female <- XB_East_West[,c(1:6)]
pi_East_male <- XB_East_West[,c(1:5,7)]
pi_West_female <- XB_East_West[,c(1:5,8)]
pi_West_male <- XB_East_West[,c(1:5,9)]
pi_Lab_female <- XB_East_West[,c(1:5,10)]
pi_Lab_male <- XB_East_West[,c(1:5,11)]
# add a faceting variable
pi_East_female$pop <- "pi_East_female"
pi_East_male$pop <- "pi_East_male"
pi_West_female$pop <- "pi_West_female"
pi_West_male$pop <- "pi_West_male"
pi_Lab_female$pop <- "pi_Lab_female"
pi_Lab_male$pop <- "pi_Lab_male"

colnames(pi_East_female) <- c("Chr","start","end","mid","sites","pi","pop")
colnames(pi_East_male) <- c("Chr","start","end","mid","sites","pi","pop")
colnames(pi_West_female) <- c("Chr","start","end","mid","sites","pi","pop")
colnames(pi_West_male) <- c("Chr","start","end","mid","sites","pi","pop")
colnames(pi_Lab_female) <- c("Chr","start","end","mid","sites","pi","pop")
colnames(pi_Lab_male) <- c("Chr","start","end","mid","sites","pi","pop")


overlay_plot_Data <- rbind(pi_East_female,pi_East_male,
                           pi_West_female,pi_West_male,
                           pi_Lab_female,pi_Lab_male)
options(scipen = 999)

png(filename = "pi_overlay_100000bp_Chr8L.png",w=1200, h=800,units = "px", bg="transparent")
p<-ggplot(overlay_plot_Data, aes(x=start/1000000, y=pi, col = pop)) + 
  #geom_point(aes(color = factor(variable))) +
  scale_color_manual(labels = c("pi_East_female","pi_East_male",
                                "pi_Lab_female","pi_Lab_male",
                                "pi_West_female","pi_West_male"
                                ),  
                     values=c("red", "pink", "blue","light blue","gray","dark gray")) +
  geom_point(size=0.25) + #, alpha = 0.2) +
  geom_line()+
  #geom_line(aes(colour = "red"), linetype = 1) +
  scale_y_continuous(name="Pairwise nucleotide diversity") + #, limits=c(-100,100)) +
  # log transform y-axis
  scale_x_continuous(name="Position on Chr8L (Mb)", breaks=seq(0,124,10)) +
  # get rid of gray background
  theme_bw() +
#theme(legend.position = "none")
  facet_wrap( ~ pop, nrow=4) +theme(strip.text.x = element_blank())
p
dev.off()

Westonly_plot_Data <- rbind(pi_West_female,pi_West_male)

##### Plot the weird section in West Female and Male #####
png(filename = "pi_West_weirdsection_100000bp.png",w=1200, h=800,units = "px", bg="transparent")
p<-ggplot(Westonly_plot_Data, aes(x=start/1000000, y=pi, col = pop)) + 
  #geom_point(aes(color = factor(variable))) +
  scale_color_manual(labels = c("pi_West_female","pi_West_male"),  
  values=c("gray","dark gray")) +
  geom_point(size=0.25) + #, alpha = 0.2) +
  geom_line()+
  #geom_line(aes(colour = "red"), linetype = 1) +
  scale_y_continuous(name="Pairwise nucleotide diversity") + #, limits=c(-100,100)) +
  # log transform y-axis
  scale_x_continuous(name="Position on Chr8L (Mb)", breaks=seq(0,124,2), limits=c(25,32)) +
  # get rid of gray background
  theme_bw() +
  #theme(legend.position = "none")
  facet_wrap( ~ pop, nrow=4) +theme(strip.text.x = element_blank())
p
dev.off()


#### DIVERGENCE ####

# make a new df of dxy for overlay plots
dxy_East_female_East_male <- XB_East_West[,c(1:5,12)]
dxy_East_female_West_female <- XB_East_West[,c(1:5,13)]
dxy_East_female_West_male <- XB_East_West[,c(1:5,14)]
dxy_East_female_Lab_female <- XB_East_West[,c(1:5,15)]
dxy_East_female_Lab_male <- XB_East_West[,c(1:5,16)]
dxy_East_male_West_female <- XB_East_West[,c(1:5,17)]
dxy_East_male_West_male <- XB_East_West[,c(1:5,18)]
dxy_East_male_Lab_female <- XB_East_West[,c(1:5,19)]
dxy_East_male_Lab_male <- XB_East_West[,c(1:5,20)]
dxy_West_female_West_male <- XB_East_West[,c(1:5,21)]
dxy_West_female_Lab_female <- XB_East_West[,c(1:5,22)]
dxy_West_female_Lab_male <- XB_East_West[,c(1:5,23)]
dxy_West_male_Lab_female <- XB_East_West[,c(1:5,24)]
dxy_West_male_Lab_male <- XB_East_West[,c(1:5,25)]
dxy_Lab_female_Lab_male <- XB_East_West[,c(1:5,26)]
# add a faceting variable
dxy_East_female_East_male$pop <- "dxy_East_female_East_male"
dxy_East_female_West_female$pop <- "dxy_East_female_West_female"
dxy_East_female_West_male$pop <- "dxy_East_female_West_male"
dxy_East_female_Lab_female$pop <- "dxy_East_female_Lab_female"
dxy_East_female_Lab_male$pop <- "dxy_East_female_Lab_male"
dxy_East_male_West_female$pop <- "dxy_East_male_West_female"
dxy_East_male_West_male$pop <- "dxy_East_male_West_male"
dxy_East_male_Lab_female$pop <- "dxy_East_male_Lab_female"
dxy_East_male_Lab_male$pop <- "dxy_East_male_Lab_male"
dxy_West_female_West_male$pop <- "dxy_West_female_West_male"
dxy_West_female_Lab_female$pop <- "dxy_West_female_Lab_female"
dxy_West_female_Lab_male$pop <- "dxy_West_female_Lab_male"
dxy_West_male_Lab_female$pop <- "dxy_West_male_Lab_female"
dxy_West_male_Lab_male$pop <- "dxy_West_male_Lab_male"
dxy_Lab_female_Lab_male$pop <- "dxy_Lab_female_Lab_male"

colnames(dxy_East_female_East_male) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_East_female_West_female) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_East_female_West_male) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_East_female_Lab_female) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_East_female_Lab_male) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_East_male_West_female) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_East_male_West_male) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_East_male_Lab_female) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_East_male_Lab_male) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_West_female_West_male) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_West_female_Lab_female) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_West_female_Lab_male) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_West_male_Lab_female) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_West_male_Lab_male) <- c("Chr","start","end","mid","sites","dxy","pop")
colnames(dxy_Lab_female_Lab_male) <- c("Chr","start","end","mid","sites","dxy","pop")


#overlay_plot_Data <- melt(temp2[,c(2,9:11)], id="Window")
overlay_plot_Data <- rbind(dxy_East_female_East_male,dxy_East_female_West_female,
                           dxy_East_female_West_male,dxy_East_female_Lab_female,
                           dxy_East_female_Lab_male,dxy_East_male_West_female,
                           dxy_East_male_West_male,dxy_East_male_Lab_female,
                           dxy_East_male_Lab_male,dxy_West_female_West_male,
                           dxy_West_female_Lab_female,dxy_West_female_Lab_male,
                           dxy_West_male_Lab_female,dxy_West_male_Lab_male,
                           dxy_Lab_female_Lab_male)

png(filename = "dxy_6pops_100000bp_Chr8L.png",w=1200, h=800,units = "px", bg="transparent")
p<-ggplot(overlay_plot_Data, aes(x=start/1000000, y=dxy, col = pop)) + 
  #geom_point(aes(color = factor(variable))) +
#  scale_color_manual(labels = c("dxy_East_female_East_male","dxy_East_female_West_female",
#                                "dxy_East_female_West_male","dxy_East_female_Lab_female",
#                                "dxy_East_female_Lab_male","dxy_East_male_West_female",
#                                "dxy_East_male_Lab_male","dxy_West_female_West_male",
#                                "dxy_West_female_Lab_female","dxy_West_female_Lab_male",
#                                "dxy_West_male_Lab_female","dxy_West_male_Lab_male",
#                                "dxy_Lab_female_Lab_male"
#                                ),  
#                     values=c("dark red","red","pink", 
#                              "dark blue","blue","light blue",
#                              "gray","dark gray","black",
#                              "orange","yellow","dark green",
#                              "light green","purple","brown")) +
  geom_point(size=0.25) + #, alpha = 0.2) +
  geom_line()+
  #geom_line(aes(colour = "red"), linetype = 1) +
  scale_y_continuous(name="Pairwise divergence") + #, limits=c(-100,100)) +
  # log transform y-axis
  scale_x_continuous(name="Position on Chr8L (Mb)", breaks=seq(0,124,10)) +
  # get rid of gray background
  theme_bw() +
  #theme(legend.position = "none")
  facet_wrap( ~ pop, nrow=6) + 
  theme(strip.text.x = element_blank())
p
dev.off()

###### FST #######

# make a new df of fst for overlay plots
fst_East_female_East_male <- XB_East_West[,c(1:5,27)]
fst_East_female_West_female <- XB_East_West[,c(1:5,28)]
fst_East_female_West_male <- XB_East_West[,c(1:5,29)]
fst_East_female_Lab_female <- XB_East_West[,c(1:5,30)]
fst_East_female_Lab_male <- XB_East_West[,c(1:5,31)]
fst_East_male_West_female <- XB_East_West[,c(1:5,32)]
fst_East_male_West_male <- XB_East_West[,c(1:5,33)]
fst_East_male_Lab_female <- XB_East_West[,c(1:5,34)]
fst_East_male_Lab_male <- XB_East_West[,c(1:5,35)]
fst_West_female_West_male <- XB_East_West[,c(1:5,36)]
fst_West_female_Lab_female <- XB_East_West[,c(1:5,37)]
fst_West_female_Lab_male <- XB_East_West[,c(1:5,38)]
fst_West_male_Lab_female <- XB_East_West[,c(1:5,39)]
fst_West_male_Lab_male <- XB_East_West[,c(1:5,40)]
fst_Lab_female_Lab_male <- XB_East_West[,c(1:5,41)]
# add a faceting variable
fst_East_female_East_male$pop <- "fst_East_female_East_male"
fst_East_female_West_female$pop <- "fst_East_female_West_female"
fst_East_female_West_male$pop <- "fst_East_female_West_male"
fst_East_female_Lab_female$pop <- "fst_East_female_Lab_female"
fst_East_female_Lab_male$pop <- "fst_East_female_Lab_male"
fst_East_male_West_female$pop <- "fst_East_male_West_female"
fst_East_male_West_male$pop <- "fst_East_male_West_male"
fst_East_male_Lab_female$pop <- "fst_East_male_Lab_female"
fst_East_male_Lab_male$pop <- "fst_East_male_Lab_male"
fst_West_female_West_male$pop <- "fst_West_female_West_male"
fst_West_female_Lab_female$pop <- "fst_West_female_Lab_female"
fst_West_female_Lab_male$pop <- "fst_West_female_Lab_male"
fst_West_male_Lab_female$pop <- "fst_West_male_Lab_female"
fst_West_male_Lab_male$pop <- "fst_West_male_Lab_male"
fst_Lab_female_Lab_male$pop <- "fst_Lab_female_Lab_male"

colnames(fst_East_female_East_male) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_East_female_West_female) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_East_female_West_male) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_East_female_Lab_female) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_East_female_Lab_male) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_East_male_West_female) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_East_male_West_male) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_East_male_Lab_female) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_East_male_Lab_male) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_West_female_West_male) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_West_female_Lab_female) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_West_female_Lab_male) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_West_male_Lab_female) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_West_male_Lab_male) <- c("Chr","start","end","mid","sites","fst","pop")
colnames(fst_Lab_female_Lab_male) <- c("Chr","start","end","mid","sites","fst","pop")


#overlay_plot_Data <- melt(temp2[,c(2,9:11)], id="Window")
overlay_plot_Data <- rbind(fst_East_female_East_male,fst_East_female_West_female,
                           fst_East_female_West_male,fst_East_female_Lab_female,
                           fst_East_female_Lab_male,fst_East_male_West_female,
                           fst_East_male_West_male,fst_East_male_Lab_female,
                           fst_East_male_Lab_male,fst_West_female_West_male,
                           fst_West_female_Lab_female,fst_West_female_Lab_male,
                           fst_West_male_Lab_female,fst_West_male_Lab_male,
                           fst_Lab_female_Lab_male)

png(filename = "fst_6pops_100000bp_Chr8L.png",w=1200, h=800,units = "px", bg="transparent")
p<-ggplot(overlay_plot_Data, aes(x=start/1000000, y=fst, col = pop)) + 
  #geom_point(aes(color = factor(variable))) +
  #  scale_color_manual(labels = c("fst_East_female_East_male","fst_East_female_West_female",
  #                                "fst_East_female_West_male","fst_East_female_Lab_female",
  #                                "fst_East_female_Lab_male","fst_East_male_West_female",
  #                                "fst_East_male_Lab_male","fst_West_female_West_male",
  #                                "fst_West_female_Lab_female","fst_West_female_Lab_male",
  #                                "fst_West_male_Lab_female","fst_West_male_Lab_male",
  #                                "fst_Lab_female_Lab_male"
  #                                ),  
  #                     values=c("dark red","red","pink", 
  #                              "dark blue","blue","light blue",
#                              "gray","dark gray","black",
#                              "orange","yellow","dark green",
#                              "light green","purple","brown")) +
geom_point(size=0.25) + #, alpha = 0.2) +
  geom_line()+
  #geom_line(aes(colour = "red"), linetype = 1) +
  scale_y_continuous(name="Pairwise divergence") + #, limits=c(-100,100)) +
  # log transform y-axis
  scale_x_continuous(name="Position on Chr8L (Mb)", breaks=seq(0,124,10)) +
  # get rid of gray background
  theme_bw() +
  #theme(legend.position = "none")
  facet_wrap( ~ pop, nrow=6) + 
  theme(strip.text.x = element_blank())
p
dev.off()













# Open a png file

png(filename = "fst_7.5mil_to_12.5mil.png",w=1200, h=800,units = "px", bg="transparent")
  p<-ggplot(overlay_plot_Data, aes(x=start, y=fst, col = comparison)) + 
    #geom_point(aes(color = factor(variable))) +
    #scale_color_manual(breaks = c("WW_minus_WZ", "WW_minus_WY", "WZ_minus_WY"), values=c("red", "blue", "gray")) +
    geom_point(size=0.25) + #, alpha = 0.2) +
    geom_line()+
    #geom_line(aes(colour = "red"), linetype = 1) +
    scale_y_continuous(name="Fst") + #, limits=c(-100,100)) +
    # log transform y-axis
    scale_x_continuous(name="Position on Chr7", limits=c(7500000, 12500000)) +
    # get rid of gray background
    theme_bw() 
    # Get rid of the legend
    #theme(legend.position = "none")
  p + facet_wrap( ~ comparison, nrow=3)
dev.off()


temp1 <- merge(XT10_WZ,XT11_WW, by = "Window")
temp2 <- merge(temp1,XT7_WY, by = "Window")

temp2$WZ_minus_WW <- temp2$pi_window.x - temp2$pi_window.y
temp2$WY_minus_WW <- temp2$pi_window - temp2$pi_window.y
temp2$WZ_minus_WY <- temp2$pi_window.x - temp2$pi_window
subset_diff <- temp2[,c(1,14:16)]
subset_diff_overlay_plot <- melt(subset_diff, id="Window")

png(filename = "difference_diversity.png",w=1200, h=800,units = "px", bg="transparent")
p<-ggplot(subset_diff_overlay_plot, aes(x=Window, y=value, col= variable)) + 
  #geom_point(aes(color = factor(variable))) +
  #scale_color_manual(breaks = c("WW_minus_WZ", "WW_minus_WY", "WZ_minus_WY"), values=c("red", "blue", "gray")) +
  geom_point(size=0.25) + #, alpha = 0.2) +
  geom_line()+
  #geom_line(aes(colour = "red"), linetype = 1) +
  scale_y_continuous(name="Pairwise nucleotide diversity") + #, limits=c(-100,100)) +
  # log transform y-axis
  scale_x_continuous(name="Position on Chr7", limits=c(9000000, 13000000)) +
  geom_hline(yintercept = 0) +
  # get rid of gray background
  theme_bw() 
# Get rid of the legend
#theme(legend.position = "none")
p + facet_wrap( ~ variable, nrow=3)
dev.off()
```
