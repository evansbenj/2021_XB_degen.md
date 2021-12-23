# WGS of wild and lab individuals
I mapped the WGS data from the F and M lab individuals from a previous study to Austin's genome.

I also have new data from 4 wild individuals - F and M from west and east

Plan:
* trimmomatic
* map
* add readgroups
* dedup
* haplotypecaller
* combineGVCFs
* GenotypeGVCFs
* VariantFiltration
* SelectVariants
* Maybe filter again for coverage

# Depth
boxplot
```R
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2021_Xborealis_sexchr_degen/WGS_depth_filter")
library (ggplot2)

dat<-read.table("./allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz.table",header=TRUE)
# subset each dataset
west_female <- dat[,c(1:3)]
west_male <- dat[,c(1,2,4)]
east_female <- dat[,c(1,2,5)]
east_male <- dat[,c(1,2,6)]
lab_male <- dat[,c(1,2,7)]
lab_female <- dat[,c(1,2,8)]
# add a sample variable
west_female$Individual <- "West female"
west_male$Individual <- "West male"
east_female$Individual <- "East female"
east_male$Individual <- "East male"
lab_female$Individual <- "Lab female"
lab_male$Individual <- "Lab male"
# fix the column names
colnames(lab_male)[3] <- "Depth"
colnames(lab_female)[3] <- "Depth"
colnames(west_male)[3] <- "Depth"
colnames(west_female)[3] <- "Depth"
colnames(east_male)[3] <- "Depth"
colnames(east_female)[3] <- "Depth"
# combine them
newdat <- rbind(east_female,east_male,lab_female,lab_male,west_female,west_male)
# make a column for SL_or_not
newdat$SL_or_not <- "light blue" # not_SL
newdat$SL_or_not[(newdat$CHROM == "Chr8L") & (newdat$POS < 54100000)] <- "red" # this is SL
# get some numbers (median depth of each section for each sample)


#east_female_SL
median(newdat$Depth[(newdat$Individual == "East female") & 
                    (newdat$CHROM == "Chr8L") 
                  & (newdat$POS < 54100000)], na.rm=T)
# [1] 48
#east_female_notSL
median(newdat$Depth[(newdat$Individual == "East female") & 
                    (newdat$CHROM == "Chr8L") 
                  & (newdat$POS > 54100000)], na.rm=T)
# [1] 46
#east_male_SL
median(newdat$Depth[(newdat$Individual == "East male") & 
                    (newdat$CHROM == "Chr8L") 
                  & (newdat$POS < 54100000)], na.rm=T)
# [1] 44
#east_male_notSL
median(newdat$Depth[(newdat$Individual == "East male") & 
                    (newdat$CHROM == "Chr8L") 
                  & (newdat$POS > 54100000)], na.rm=T)
# [1] 43
#west_female_SL
median(newdat$Depth[(newdat$Individual == "West female") & 
                    (newdat$CHROM == "Chr8L") 
                  & (newdat$POS < 54100000)], na.rm=T)
# [1] 32
#west_female_notSL
median(newdat$Depth[(newdat$Individual == "West female") & 
                    (newdat$CHROM == "Chr8L") 
                  & (newdat$POS > 54100000)], na.rm=T)
# [1] 31

# make the SL first
newdat$SL_or_not <- factor(newdat$SL_or_not,
                           levels = c("red","light blue"),ordered = TRUE)

png("SL_notSL_depth_by_individual_Chr8L.png", width = 200, height = 110, units='mm', res = 100) 
ggplot(newdat, aes(x=Individual, y=Depth, fill=SL_or_not)) + 
  #geom_hline(yintercept=0, linetype="solid",color = "black", size=0.5) + 
  geom_boxplot(outlier.size = 0.5) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_classic() +
  scale_fill_manual(values = c("red", "light blue"),
                    name="<54.1 Mb?",
                    breaks=c("red", "light blue"),
                    labels=c("Yes", "No")) +
  theme(axis.text=element_text(size=12),
         axis.title=element_text(size=14))+
  scale_y_continuous(name="Depth (Chr8L)")+
  theme(legend.position = c(0.95, 0.8))
dev.off()




# Chr7S
# clear the environment
rm(list = ls())
dat<-read.table("./allsites_Chr7S_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz.table",header=TRUE)
# subset each dataset
west_female <- dat[,c(1:3)]
west_male <- dat[,c(1,2,4)]
east_female <- dat[,c(1,2,5)]
east_male <- dat[,c(1,2,6)]
lab_male <- dat[,c(1,2,7)]
lab_female <- dat[,c(1,2,8)]
# add a sample variable
west_female$Individual <- "West female"
west_male$Individual <- "West male"
east_female$Individual <- "East female"
east_male$Individual <- "East male"
lab_female$Individual <- "Lab female"
lab_male$Individual <- "Lab male"
# fix the column names
colnames(lab_male)[3] <- "Depth"
colnames(lab_female)[3] <- "Depth"
colnames(west_male)[3] <- "Depth"
colnames(west_female)[3] <- "Depth"
colnames(east_male)[3] <- "Depth"
colnames(east_female)[3] <- "Depth"
# combine them
newdat <- rbind(east_female,east_male,lab_female,lab_male,west_female,west_male)
# make a column for SL_or_not
newdat$SL_or_not <- "light blue" # not_SL
newdat$SL_or_not[(newdat$CHROM == "Chr7S") & (newdat$POS < 23000000)] <- "red" # this is SL
# get some numbers (median depth of each section for each sample)


#east_female_SL
median(newdat$Depth[(newdat$Individual == "East female") & 
                      (newdat$CHROM == "Chr7S") 
                    & (newdat$POS < 23000000)], na.rm=T)
# [1] 45 (48 for chr8L)
#east_female_notSL
median(newdat$Depth[(newdat$Individual == "East female") & 
                      (newdat$CHROM == "Chr7S") 
                    & (newdat$POS > 23000000)], na.rm=T)
# [1] 45 (46 for chr8L)
#east_male_SL
median(newdat$Depth[(newdat$Individual == "East male") & 
                      (newdat$CHROM == "Chr7S") 
                    & (newdat$POS < 23000000)], na.rm=T)
# [1] 43 (44 for chr8L)
#east_male_notSL
median(newdat$Depth[(newdat$Individual == "East male") & 
                      (newdat$CHROM == "Chr7S") 
                    & (newdat$POS > 23000000)], na.rm=T)
# [1] 43 (43 for chr8L)
#west_female_SL
median(newdat$Depth[(newdat$Individual == "West female") & 
                      (newdat$CHROM == "Chr7S") 
                    & (newdat$POS < 23000000)], na.rm=T)
# [1] 30 (32 for chr8L)
#west_female_notSL
median(newdat$Depth[(newdat$Individual == "West female") & 
                      (newdat$CHROM == "Chr7S") 
                    & (newdat$POS > 23000000)], na.rm=T)
# [1] 31 (31 for chr8L)

# make the SL first
newdat$SL_or_not <- factor(newdat$SL_or_not,
                           levels = c("red","light blue"),ordered = TRUE)

png("SL_notSL_depth_by_individual_Chr7S.png", width = 200, height = 110, units='mm', res = 100) 
ggplot(newdat, aes(x=Individual, y=Depth, fill=SL_or_not)) + 
  geom_boxplot(outlier.size = 0.5) +
  coord_cartesian(ylim = c(0, 100)) +
  theme_classic() +
  scale_fill_manual(values = c("red", "light blue"),
                    name="<23.0 Mb?",
                    breaks=c("red", "light blue"),
                    labels=c("Yes", "No")) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  scale_y_continuous(name="Depth (Chr7S)")
dev.off()
```

```R
library (ggplot2)
library(reshape2) # this facilitates the overlay plot
library(tidyverse)
library(tidyquant) # needed for moving averages
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2021_Xborealis_sexchr_degen/WGS_depth_filter")
# get all the files with the site data
depth <- read.table("allsites_Chr7S_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz.table", header = T, sep="\t")
overlay_plot_Data <- melt(depth, id=c("CHROM","POS"))
colnames(overlay_plot_Data)[3] <- "individual"
colnames(overlay_plot_Data)[4] <- "Depth"
# Make a new column with nicer individual names
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
        "BJE4441_female_west_merge_sorted_dedup.bam.DP")] <- "West female" 
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
        "BJE4442_male_west_merge_sorted_dedup.bam.DP")] <- "West male" 
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
        "BJE4515_female_east_merge_sorted_dedup.bam.DP")] <- "East female" 
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
        "BJE4536_male_east_merge_sorted_dedup.bam.DP")] <- "East male" 
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
        "SRR6357672.1_trim_fixed_merged_dedup.bam.DP")] <- "Lab male" 
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
        "SRR6357673_trim_sorted_dedup.bam.DP")] <- "Lab female" 

# data are now read for overlay plotting
dim(overlay_plot_Data)
dim(depth)

png(filename = "Depth_bysite_6samples_Chr7S.png",w=1200, h=800,units = "px", bg="transparent")
p<-ggplot(overlay_plot_Data, aes(x=POS/1000000, y=Depth, col = Individual)) + 
  geom_ma(ma_fun = SMA, n = 1000) +
  #geom_line(aes(colour = "red"), linetype = 1) +
  scale_y_continuous(name="Depth", limits=c(0,49)) +
  scale_x_continuous(name="Position on Chr7S (Mb)", breaks=seq(0,124,10)) +
  # get rid of gray background
  theme_bw() +
  #theme(legend.position = "none")
  facet_wrap( ~ Individual, nrow=6) + 
  theme(strip.text.x = element_blank())+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20))+ 
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=16))
p
dev.off()

# clear the environment
rm(list = ls())
# get all the files with the site data
depth <- read.table("allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz.table", header = T, sep="\t")
overlay_plot_Data <- melt(depth, id=c("CHROM","POS"))
colnames(overlay_plot_Data)[3] <- "individual"
colnames(overlay_plot_Data)[4] <- "Depth"
# Make a new column with nicer individual names
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
                                "BJE4441_female_west_merge_sorted_dedup.bam.DP")] <- "west female" 
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
                                "BJE4442_male_west_merge_sorted_dedup.bam.DP")] <- "west male" 
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
                                "BJE4515_female_east_merge_sorted_dedup.bam.DP")] <- "east female" 
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
                                "BJE4536_male_east_merge_sorted_dedup.bam.DP")] <- "east male" 
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
                                "SRR6357672.1_trim_fixed_merged_dedup.bam.DP")] <- "lab male" 
overlay_plot_Data$Individual[(overlay_plot_Data$individual == 
                                "SRR6357673_trim_sorted_dedup.bam.DP")] <- "lab female" 

# data are now read for overlay plotting
dim(overlay_plot_Data)
dim(depth)

png(filename = "Depth_bysite_6samples_Chr8L.png",w=1200, h=800,units = "px", bg="transparent")
p<-ggplot(overlay_plot_Data, aes(x=POS/1000000, y=Depth, col = Individual)) + 
  geom_ma(ma_fun = SMA, n = 1000) +
  #geom_line(aes(colour = "red"), linetype = 1) +
  scale_y_continuous(name="Depth", limits=c(0,49)) +
  scale_x_continuous(name="Position on Chr8L (Mb)", breaks=seq(0,124,10)) +
  # get rid of gray background
  theme_bw() +
  #theme(legend.position = "none")
  facet_wrap( ~ Individual, nrow=6) + 
  theme(strip.text.x = element_blank())+
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20))+ 
  theme(legend.title=element_text(size=20), 
        legend.text=element_text(size=16))
p
dev.off()
```
