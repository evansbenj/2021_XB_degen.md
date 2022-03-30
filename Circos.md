# Circos plot of sex linkage
```R
## Working directory
setwd("/Users/Shared/Previously Relocated Items/Security/projects/submitted/2021_Xborealis_sexchr_degen/Circos_plot")
library(tidyverse)
library(ggplot2)
library(circlize)
library(dplyr)
library(stringr)
# https://jokergoo.github.io/circlize_book/book/initialize-genomic-plot.html
# https://cran.r-project.org/web/packages/circlize/circlize.pdf
# initilize
rm(list=ls()) # removes all variables

sample_vector <- c("east_only","XBlab","west_only")
chrs <- factor(c("Chr1L","Chr2L","Chr3L","Chr4L","Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                 "Chr1S","Chr2S","Chr3S","Chr4S","Chr5S","Chr6S","Chr7S","Chr8S","Chr9_10S"))

#png(paste("XB_circos",eval(analysis),".png",sep=""),
#    width = 500, height = 500, units='mm', res = 100)
# https://jokergoo.github.io/circlize_book/book/introduction.html
# Initialize library

png("sex_linkage_circos.png",
    width = 100, height = 100, units='mm', res = 300)

# set up chr structure
chr_structure <- data.frame(CHR= c("Chr1L","Chr2L","Chr3L","Chr4L","Chr5L",
                                   "Chr6L","Chr7L","Chr8L","Chr9_10L",
                                   "Chr9_10S","Chr8S","Chr7S","Chr6S","Chr5S",
                                   "Chr4S","Chr3S","Chr2S","Chr1S"),
                           start = c(rep(1,18)),
                           end = c(232529968,184566230,145564450,
                                   156120766,174499025,157843503,
                                   136892545,123836260,135078615,
                                   110702965,105436523,105895007,
                                   137668414,139053355,131359389,
                                   127416163,167897112,196169797))
# Make the chr an ordered factor - this is necessary to ensure
# proper initialization of circos plot
chr_structure$CHR <- factor(chr_structure$CHR,
                            levels = c("Chr1L","Chr2L","Chr3L","Chr4L","Chr5L",
                                   "Chr6L","Chr7L","Chr8L","Chr9_10L",
                                   "Chr9_10S","Chr8S","Chr7S","Chr6S","Chr5S",
                                   "Chr4S","Chr3S","Chr2S","Chr1S"),
                              ordered = TRUE)
# initialize plot with custom genome
# https://mran.microsoft.com/snapshot/2014-12-11/web/packages/circlize/vignettes/genomic_plot.pdf
circos.clear()
circos.par("gap.degree" = c(rep(3, 18)), 
           "cell.padding" = c(0, 0, 0, 0), 
           "start.degree" = 90,
           "track.height" = 0.15) 

circos.initializeWithIdeogram(chr_structure, plotType = c("axis", "labels"),
                              sort.chr = F,
                              chromosome.index = c("Chr1L","Chr2L","Chr3L","Chr4L","Chr5L",
                                                   "Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr9_10S","Chr8S","Chr7S","Chr6S",
                                                   "Chr5S","Chr4S","Chr3S","Chr2S","Chr1S"),
                              tickLabelsStartFromZero = FALSE,
                              major.by = 100000000,
                              axis.labels.cex = 0.3,
                              labels.cex = 0.5)


#chr8l
draw.sector(286, 293, rou1 = 0.87, rou2 = 0.57, 
            clock.wise = FALSE, col = "lightpink", border = NA) # red over track 1

#chr7S
draw.sector(223, 226.5, rou1 = 0.87, rou2 = 0.69, 
            clock.wise = FALSE, col = "cadetblue1", border = NA) # red over track 1

# center L
draw.sector(259, 91, rou1 = 0.35, rou2 = 0.0, 
            clock.wise = FALSE, col = "whitesmoke", border = NA) # red over track 1

# loop through each sample
for (sample in sample_vector){
  a <- read.table(paste(eval(sample), "_filtered.assoc", sep=""), header=TRUE)
  #a <- read.table("./east_only_filtered.assoc",header=TRUE)
  # make the chr a factor so we can use it for faceting
  a$CHR <- as.factor(a$CHR)
  # make a new column which is -logP
  a$minus_log_P<- -log(a$P)
  # make a dataframe for circular plotting
  Allchr_circular <- as.data.frame(a[,c(1,3,9,11)])
  # make a column to color values with low pvalues
  Allchr_circular$color <- "light gray"
  Allchr_circular[c("color")][which(Allchr_circular$P < 0.003), ] <- 'red'
  Allchr_circular[c("color")][which((Allchr_circular$P > 0.003)&(Allchr_circular$P < 0.03)), ] <- 'orange'
  
  # generate a variable that determines which points are plotted first
  Allchr_circular$order <- 1
  Allchr_circular[c("order")][which(Allchr_circular$P < 0.003), ] <- 3
  Allchr_circular[c("order")][which((Allchr_circular$P > 0.003)&(Allchr_circular$P < 0.03)), ] <- 2
  
  # make the BP and minusP columns numeric
  Allchr_circular$BP <- as.numeric(Allchr_circular$BP)
  Allchr_circular$minus_log_P <- as.numeric(Allchr_circular$minus_log_P)
  # add an end column, rename BP column, and change order of colums:
  Allchr_circular$end <- Allchr_circular$BP
  colnames(Allchr_circular)[1] <- "chr"
  colnames(Allchr_circular)[2] <- "start"
  Allchr_circular <- Allchr_circular[,c(1,2,7,4,5)]
  # order the columns
  Allchr_circular$chr <- factor(Allchr_circular$chr, 
                                   levels = c("Chr1L","Chr2L","Chr3L","Chr4L","Chr5L","Chr6L","Chr7L","Chr8L","Chr9_10L",
                                                   "Chr9_10S","Chr8S","Chr7S","Chr6S","Chr5S","Chr4S","Chr3S","Chr2S","Chr1S"),
                                   ordered = TRUE)
  # sort the df so that the chrs are in the correct order
  Allchr_circular <- Allchr_circular[order(Allchr_circular$chr),]
  
  if (sample == 'XBlab'){
    circos.track(ylim = c(0, 50))
  } else if (sample == 'east_only'){
    circos.track(ylim = c(0, 10))
  } else {
    circos.track(ylim = c(0, 15))
  }
  circos.trackPoints(Allchr_circular$chr, x = Allchr_circular$start, 
                     y = Allchr_circular$minus_log_P, 
                     track.index = get.cell.meta.data("track.index"),
                     pch = 16, cex=0.4, col = Allchr_circular$color)

}

# label the species in the center
  text(0, 0, "S   L", cex = 1.5)

dev.off()
```
