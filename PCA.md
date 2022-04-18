# Chr8L WGS data
Working directory:
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/WGS_vcfs_by_chr/combined_andfiltereds_gvcfs/
```
Subset SL and not Sl regions
SL first:
```
bcftools view allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz --regions Chr8L:1-54100000 > allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_SLonly.vcf.gz
```
Not AL:
```
bcftools view allsites_allchrs_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz -R ../../XB_RADseq_mapped_to_AustinXB_BJE_2021/combined_and_genotyped_vcfs_trim/big_chrs_notSL.bed -O z -o allsites_allchrs_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_noSL_noChr8L_1_54100000.vcf.gz
```

R script:
```
# Note that a bug in the SNPrelate package means you need to quit R each time you reload
# a vcf file


#https://github.com/zhengxwen/SNPRelate/issues/13
#http://corearray.sourceforge.net/tutorials/SNPRelate/
#https://github.com/zhengxwen/SNPRelate/wiki/Preparing-Data

#GenotypeVCFs_noBSQR_filtered_aDNA_only.vcf.gz
#GenotypeVCFs_noBSQR_filtered_xDNA_only.vcf.gz
#GenotypeVCFs_noBSQR_filtered_aDNA_only_no_lowcoverage_individuals.vcf
setwd('/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/WGS_vcfs_by_chr/combined_andfiltereds_gvcfs')
library("devtools")
#    if (!require("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#    BiocManager::install("gdsfmt")
library(gdsfmt)
#    if (!require("BiocManager", quietly = TRUE))
#        install.packages("BiocManager")
#    BiocManager::install("SNPRelate")
library(SNPRelate)

#vcf.fn <- "XBwild_unfiltered_allchrs.vcf.gz_filtered.vcf.gz_filtered_removed_Chr7S_notSL.recode.vcf"
#vcf.fn <- "XBwild_unfiltered_allchrs.vcf.gz_filtered.vcf.gz_filtered_removed_Chr7S_XMXVXB_notSL.recode.vcf"
#vcf.fn <- "allsites_Chr7S_XMXVXB_notSL.recode.vcf"
pseudoonly <- "allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_pseudoonly.vcf.gz"
#vcf.fn <- "allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_30000000_40000000.vcf.gz"
#vcf.fn <- "allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_27000000_30000000.vcf.gz"
#vcf.fn <- "allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_63000000_93000000.vcf.gz"
SLonly <-  "allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_SLonly.vcf.gz"
#vcf.fn <- "allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed_1_27000000.vcf.gz"
snpgdsVCF2GDS(pseudoonly, "pseudo", method="biallelic.only",ignore.chr.prefix = c("chr","Scaffold","chr9_"))
snpgdsVCF2GDS(SLonly, "SL", method="biallelic.only",ignore.chr.prefix = c("chr","Scaffold","chr9_"))

pseudo_genofile = snpgdsOpen("pseudo", readonly=FALSE)
SL_genofile = snpgdsOpen("SL", readonly=FALSE)

samp.annot<-data.frame(pop.group = c("female_W","male_W","female_E","male_E","male_lab",
                                     "female_lab"))

add.gdsn(pseudo_genofile, "sample.annot", val=samp.annot)
add.gdsn(SL_genofile, "sample.annot", val=samp.annot)

snpgdsSummary("pseudo")
snpgdsSummary("SL")
#pca <- snpgdsPCA(genofile, num.thread=2)
#pc.percent <- pca$varprop*100
#head(round(pc.percent, 2))

# LD prunning I am setting ld.threshold=1 to minimize ld filtering
pseudo_snpset <- snpgdsLDpruning(pseudo_genofile, ld.threshold=0.1,  method = c("composite"), missing.rate=0, ver
bose = TRUE)
# 335 SNPs retained
SL_snpset <- snpgdsLDpruning(SL_genofile, ld.threshold=0.1,  method = c("composite"), missing.rate=0, verbose = T
RUE)
# 311 SNPs retained
pseudo_snpset.id <- unlist(pseudo_snpset)
SL_snpset.id <- unlist(SL_snpset)

# for LD prunning, uncomment next line and comment the one after that
pseudo_pca <- snpgdsPCA(pseudo_genofile, snp.id=pseudo_snpset.id, num.thread=2)
SL_pca <- snpgdsPCA(SL_genofile, snp.id=SL_snpset.id, num.thread=2)
# no LD pruning
#pca <- snpgdsPCA(genofile)
pseudo_pc.percent <- pseudo_pca$varprop*100
SL_pc.percent <- SL_pca$varprop*100
head(round(pseudo_pc.percent, 2))
# updated for WGS
# [1] 37.49 27.34 22.54 10.19  2.45  0.00
head(round(SL_pc.percent, 2))
# updated for WGS
# [1] 36.07 25.28 19.74 14.73  4.19  0.00


tab <- data.frame(sample.id = pseudo_pca$sample.id,
    EV1 = pseudo_pca$eigenvect[,1],    # the first eigenvector
    EV2 = pseudo_pca$eigenvect[,2],    # the second eigenvector
    EV3 = pseudo_pca$eigenvect[,3],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

library(ggplot2)
library(ggrepel)

tab$Species <- c("female","male","female","male","male","female")
tab$samp.color <- c("blue","blue","red","red","gray","gray")
tab$samp.fieldid <- c("female","male","female","male","male","female")
pseudo_plot<-ggplot(data=tab, aes(x=EV1,y=EV2, label = samp.fieldid, color = samp.color)) +
    # label axis 
    labs(x=expression("Eigenvector 1"), y=element_blank()) +
    # legend details
    scale_colour_manual(name="Population", values = c(
                                                    "blue"="blue",
                                                    "gray"="gray", 
                                                       "red"="red"),
                        breaks=c("blue","gray","red"),
                        labels=c("West","Lab","East"))+
    # add points and fieldID labels
    geom_text_repel(aes(EV1,EV2, label=(samp.fieldid))) + geom_point(size=4) + 
    # change to cleaner theme
    theme_classic(base_size = 16) +
    # make it clean
    theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
    # italicize species names
    #theme(legend.text = element_text(face="italic"))+ 
    # move the legend
   # theme(legend.position = c(.8, .1)) +
    # add a title
    ggtitle("Not sex-linked") + 
    # remove boxes around legend symbols
    theme(legend.key = element_blank())
pseudo_plot

SL_tab <- data.frame(sample.id = SL_pca$sample.id,
                  EV1 = SL_pca$eigenvect[,1],    # the first eigenvector
                  EV2 = SL_pca$eigenvect[,2],    # the second eigenvector
                  EV3 = SL_pca$eigenvect[,3],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

SL_tab$Species <- c("female","male","female","male","male","female")
SL_tab$samp.color <- c("blue","blue","red","red","gray","gray")
SL_tab$samp.fieldid <- c("female","male","female","male","male","female")
SL_plot<-ggplot(data=SL_tab, aes(x=EV1,y=EV2, label = samp.fieldid, color = samp.color)) +
    # label axis 
    labs(x=expression("Eigenvector 1"), y=expression("Eigenvector 2")) +
    # legend details
    scale_colour_manual(name="Population", values = c(
        "blue"="blue",
        "gray"="gray",
        "red"="red"),
        breaks=c("blue","gray","red"),
        labels=c("West","Lab","East"))+
    # add points and fieldID labels
    geom_text_repel(aes(EV1,EV2, label=(samp.fieldid))) + geom_point(size=4) + 
    # change to cleaner theme
    theme_classic(base_size = 16) +
    # make it clean
    theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
    # italicize species names
    #theme(legend.text = element_text(face="italic"))+ 
    # move the legend
    # theme(legend.position = c(.8, .1)) +
    # add a title
    ggtitle("Sex-linked (Chr8L <54.1 Mb)") + 
    # remove boxes around legend symbols
    theme(legend.position = "none")
    #theme(legend.key = element_blank())
SL_plot

#library(ggpubr)
library(gridExtra)
pdf("PCA_WGS_data_SL_and_nonSL_E1E2.pdf",w=8, h=3.6, version="1.4", bg="transparent")
      grid.arrange(SL_plot, pseudo_plot, nrow = 1, widths = c(1, 1.2))
dev.off()
```

# Chr7S
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
# Note that a bug in the SNPrelate package means you need to quit R each time you reload
# a vcf file


#https://github.com/zhengxwen/SNPRelate/issues/13
#http://corearray.sourceforge.net/tutorials/SNPRelate/
#https://github.com/zhengxwen/SNPRelate/wiki/Preparing-Data

#GenotypeVCFs_noBSQR_filtered_aDNA_only.vcf.gz
#GenotypeVCFs_noBSQR_filtered_xDNA_only.vcf.gz
#GenotypeVCFs_noBSQR_filtered_aDNA_only_no_lowcoverage_individuals.vcf
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2021_Xborealis_sexchr_degen/PCA')
library("devtools")
library(gdsfmt)
library(SNPRelate)



vcf.fn <- "XBwild_unfiltered_allchrs.vcf.gz_filtered.vcf.gz_filtered_removed_Chr7S_notSL.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only",ignore.chr.prefix = c("chr","Scaffold","chr9_"))

genofile = snpgdsOpen("test.gds", readonly=FALSE)
samp.annot<-data.frame(pop.group = c("W","W","W","W","W","W",
                                     "W","W","W","W","W","W",
                                     "W","W","W","W","C","C",
                                     "C","C","C","C","C","C",
                                     "E","E","E","E","E","W",
                                     "W","W","W","W","W","W",
                                     "W","W","W","W","W","W",
                                     "W","W","C","C","C","C",
                                     "C","E","E","E","E","E"))


add.gdsn(genofile, "sample.annot", val=samp.annot)

snpgdsSummary("test.gds")

#pca <- snpgdsPCA(genofile, num.thread=2)
#pc.percent <- pca$varprop*100
#head(round(pc.percent, 2))

# LD prinnung
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,  method = c("composite"),missing.rate=0, verbose = TRUE)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))




tab <- data.frame(sample.id = pca$sample.id,
    EV1 = pca$eigenvect[,1],    # the first eigenvector
    EV2 = pca$eigenvect[,2],    # the second eigenvector
    EV3 = pca$eigenvect[,3],    # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

#plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
#text(tab$EV2, tab$EV1,labels=tab$sample.id, cex= 0.4)

library(ggplot2)
#ggplot(...)+...+ theme(axis.text.x = element_text(angle=60, hjust=1))
#devtools::install_github("slowkow/ggrepel")
library(ggrepel)

pdf("PCA_plot_XB_W_E_Chr7S_gt23Mb.pdf",w=8, h=8, version="1.4", bg="transparent")
tab$Species <- c("W","W","W","W","W","W",
                 "W","W","W","W","W","W",
                 "W","W","W","W","C","C",
                 "C","C","C","C","C","C",
                 "E","E","E","E","E","W",
                 "W","W","W","W","W","W",
                 "W","W","W","W","W","W",
                 "W","W","C","C","C","C",
                 "C","E","E","E","E","E")
tab$samp.color <- c("blue","blue","blue","blue","blue","blue","blue","blue",
                    "blue","blue","blue","blue","blue","blue","blue","blue",
                    "gray","gray","gray","gray","gray","gray","gray","gray",
                    "red","red","red","red","red","blue","blue","blue","blue",
                    "blue","blue","blue","blue","blue","blue","blue","blue",
                    "blue","blue","blue","gray","gray","gray","gray","gray",
                    "red","red","red","red","red")
tab$samp.fieldid <- c("Fem_Cheplaskei_BJE4459","Fem_Cheplaskei_BJE4461",
                      "Fem_Cheplaskei_BJE4470","Fem_Chesuwe_BJE4479",
                      "Fem_Chesuwe_BJE4481","Fem_Eldoret_BJE4471",
                      "Fem_Eldoret_BJE4472","Fem_Eldoret_BJE4474",
                      "Fem_Eldoret_BJE4475","Fem_Eldoret_BJE4476",
                      "Fem_Kiminini_BJE4429","Fem_Kiminini_BJE4433",
                      "Fem_Lukhome_BJE4441","Fem_Lukhome_BJE4444",
                      "Fem_Lukhome_BJE4445","Fem_Lukhome_BJE4446",
                      "Fem_NMK_BJE4562","Fem_NMK_BJE4564","Fem_NMK_BJE4567",
                      "Fem_NMK_BJE4568","Fem_Nakuru_BJE4363","Fem_Nakuru_BJE4364",
                      "Fem_Nakuru_BJE4367","Fem_Nakuru_BJE4368","Fem_Wundanyi_BJE4515",
                      "Fem_Wundanyi_BJE4516","Fem_Wundanyi_BJE4534",
                      "Fem_Wundanyi_BJE4535","Fem_Wundanyi_BJE4541",
                      "Mal_Cheplaskei_BJE4460","Mal_Cheplaskei_BJE4462",
                      "Mal_Cheplaskei_BJE4465","Mal_Cheplaskei_BJE4469",
                      "Mal_Chesuwe_BJE4477","Mal_Chesuwe_BJE4478",
                      "Mal_Chesuwe_BJE4480","Mal_Eldoret_BJE4473","Mal_Kiminini_BJE4430",
                      "Mal_Kiminini_BJE4431","Mal_Kiminini_BJE4432",
                      "Mal_Kisumu_BJE4391","Mal_Lukhome_BJE4442","Mal_Lukhome_BJE4443",
                      "Mal_Lukhome_BJE4447","Mal_NMK_BJE4563","Mal_Nakuru_BJE4365",
                      "Mal_Nakuru_BJE4366","Mal_Nakuru_BJE4374","Mal_Nakuru_BJE4377",
                      "Mal_Wundanyi_BJE4536","Mal_Wundanyi_BJE4537",
                      "Mal_Wundanyi_BJE4538","Mal_Wundanyi_BJE4539",
                      "Mal_Wundanyi_BJE4540")
# red is East
# gray is Central
# blue is West
d<-ggplot(data=tab, aes(x=EV1,y=-EV2, label = samp.fieldid, color = samp.color)) +
    # label axis 
    labs(x=expression("Eigenvector 1"), y=expression("Eigenvector 2")) +
    # legend details
    scale_colour_manual(name="Population", values = c("gray"="gray", 
                                                              "blue"="blue","red"="red"),
                        breaks=c("gray","blue","red"),labels=c("Central", "West", "East"))+
    # add points and fieldID labels
    geom_text_repel(aes(EV1,-EV2, label=(samp.fieldid))) + geom_point(size=4) + 
    # change to cleaner theme
    theme_classic(base_size = 16) +
    # make it clean
    theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
    # italicize species names
    #theme(legend.text = element_text(face="italic"))+ 
    # move the legend
    theme(legend.position = c(.8, .1)) +
    # add a title
    ggtitle("Principal Components Analsis") + 
    # remove boxes around legend symbols
    theme(legend.key = element_blank())
d
dev.off()









#https://github.com/zhengxwen/SNPRelate/issues/13
#http://corearray.sourceforge.net/tutorials/SNPRelate/
#https://github.com/zhengxwen/SNPRelate/wiki/Preparing-Data

#GenotypeVCFs_noBSQR_filtered_aDNA_only.vcf.gz
#GenotypeVCFs_noBSQR_filtered_xDNA_only.vcf.gz
#GenotypeVCFs_noBSQR_filtered_aDNA_only_no_lowcoverage_individuals.vcf
setwd('/Users/Shared/Previously\ Relocated\ Items/Security/projects/2021_Xborealis_sexchr_degen/PCA')
library("devtools")
library(gdsfmt)
library(SNPRelate)



vcf.fn <- "XBwild_unfiltered_allchrs.vcf.gz_filtered.vcf.gz_filtered_removed_Chr7S_SL.recode.vcf"
snpgdsVCF2GDS(vcf.fn, "test.gds", method="biallelic.only",ignore.chr.prefix = c("chr","Scaffold","chr9_"))

genofile = snpgdsOpen("test.gds", readonly=FALSE)
samp.annot<-data.frame(pop.group = c("W","W","W","W","W","W",
                                     "W","W","W","W","W","W",
                                     "W","W","W","W","C","C",
                                     "C","C","C","C","C","C",
                                     "E","E","E","E","E","W",
                                     "W","W","W","W","W","W",
                                     "W","W","W","W","W","W",
                                     "W","W","C","C","C","C",
                                     "C","E","E","E","E","E"))


add.gdsn(genofile, "sample.annot", val=samp.annot)

snpgdsSummary("test.gds")

#pca <- snpgdsPCA(genofile, num.thread=2)
#pc.percent <- pca$varprop*100
#head(round(pc.percent, 2))

# LD pruning
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2,  method = c("composite"),missing.rate=0, verbose = TRUE)
snpset.id <- unlist(snpset)
pca <- snpgdsPCA(genofile, snp.id=snpset.id, num.thread=2)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))




tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  EV3 = pca$eigenvect[,3],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)

#plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")
#text(tab$EV2, tab$EV1,labels=tab$sample.id, cex= 0.4)

library(ggplot2)
#ggplot(...)+...+ theme(axis.text.x = element_text(angle=60, hjust=1))
#devtools::install_github("slowkow/ggrepel")
library(ggrepel)

pdf("PCA_plot_XB_W_E_Chr7S_lt23Mb.pdf",w=8, h=8, version="1.4", bg="transparent")
tab$Species <- c("W","W","W","W","W","W",
                 "W","W","W","W","W","W",
                 "W","W","W","W","C","C",
                 "C","C","C","C","C","C",
                 "E","E","E","E","E","W",
                 "W","W","W","W","W","W",
                 "W","W","W","W","W","W",
                 "W","W","C","C","C","C",
                 "C","E","E","E","E","E")
tab$samp.color <- c("blue","blue","blue","blue","blue","blue","blue","blue",
                    "blue","blue","blue","blue","blue","blue","blue","blue",
                    "gray","gray","gray","gray","gray","gray","gray","gray",
                    "red","red","red","red","red","blue","blue","blue","blue",
                    "blue","blue","blue","blue","blue","blue","blue","blue",
                    "blue","blue","blue","gray","gray","gray","gray","gray",
                    "red","red","red","red","red")
tab$samp.fieldid <- c("Fem_Cheplaskei_BJE4459","Fem_Cheplaskei_BJE4461",
                      "Fem_Cheplaskei_BJE4470","Fem_Chesuwe_BJE4479",
                      "Fem_Chesuwe_BJE4481","Fem_Eldoret_BJE4471",
                      "Fem_Eldoret_BJE4472","Fem_Eldoret_BJE4474",
                      "Fem_Eldoret_BJE4475","Fem_Eldoret_BJE4476",
                      "Fem_Kiminini_BJE4429","Fem_Kiminini_BJE4433",
                      "Fem_Lukhome_BJE4441","Fem_Lukhome_BJE4444",
                      "Fem_Lukhome_BJE4445","Fem_Lukhome_BJE4446",
                      "Fem_NMK_BJE4562","Fem_NMK_BJE4564","Fem_NMK_BJE4567",
                      "Fem_NMK_BJE4568","Fem_Nakuru_BJE4363","Fem_Nakuru_BJE4364",
                      "Fem_Nakuru_BJE4367","Fem_Nakuru_BJE4368","Fem_Wundanyi_BJE4515",
                      "Fem_Wundanyi_BJE4516","Fem_Wundanyi_BJE4534",
                      "Fem_Wundanyi_BJE4535","Fem_Wundanyi_BJE4541",
                      "Mal_Cheplaskei_BJE4460","Mal_Cheplaskei_BJE4462",
                      "Mal_Cheplaskei_BJE4465","Mal_Cheplaskei_BJE4469",
                      "Mal_Chesuwe_BJE4477","Mal_Chesuwe_BJE4478",
                      "Mal_Chesuwe_BJE4480","Mal_Eldoret_BJE4473","Mal_Kiminini_BJE4430",
                      "Mal_Kiminini_BJE4431","Mal_Kiminini_BJE4432",
                      "Mal_Kisumu_BJE4391","Mal_Lukhome_BJE4442","Mal_Lukhome_BJE4443",
                      "Mal_Lukhome_BJE4447","Mal_NMK_BJE4563","Mal_Nakuru_BJE4365",
                      "Mal_Nakuru_BJE4366","Mal_Nakuru_BJE4374","Mal_Nakuru_BJE4377",
                      "Mal_Wundanyi_BJE4536","Mal_Wundanyi_BJE4537",
                      "Mal_Wundanyi_BJE4538","Mal_Wundanyi_BJE4539",
                      "Mal_Wundanyi_BJE4540")
# red is East
# gray is Central
# blue is West
d<-ggplot(data=tab, aes(x=EV1,y=-EV2, label = samp.fieldid, color = samp.color)) +
    # label axis 
    labs(x=expression("Eigenvector 1"), y=expression("Eigenvector 2")) +
    # legend details
    scale_colour_manual(name="Species/Population", values = c("gray"="gray", "blue"="blue","red"="red"),breaks=c("gray","blue","red"),labels=c("X. victorianus", "X. borealis West", "X. borealis East"))+
    # add points and fieldID labels
    geom_text_repel(aes(EV1,-EV2, label=(samp.fieldid))) + geom_point(size=4) + 
    # change to cleaner theme
    theme_classic(base_size = 16) +
    # make it clean
    theme_bw()+ theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank()) + 
    # italicize species names
    theme(legend.text = element_text(face="italic"))+ 
    # move the legend
    theme(legend.position = c(.1, .1)) +
    # add a title
    ggtitle("Principal Components Analsis") + 
    # remove boxes around legend symbols
    theme(legend.key = element_blank())
d
dev.off()

```
