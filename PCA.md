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
