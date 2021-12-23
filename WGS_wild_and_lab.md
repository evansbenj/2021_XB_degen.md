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
