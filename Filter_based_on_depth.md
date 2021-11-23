# Filtering based on rolling average of depth

Get depth
```
#!/bin/sh
#SBATCH --job-name=VariantsToTable
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=10gb
#SBATCH --output=VariantToTable.%J.out
#SBATCH --error=VariantToTable.%J.err
#SBATCH --account=def-ben


# This script will execute the GATK command "VariantsToTable

# execute like this:
# sbatch 2021_VariantFiltration.sh path_and_file

module load nixpkgs/16.09 gatk/4.1.0.0

gatk --java-options -Xmx8G VariantsToTable -V ${1} -F CHROM -F POS -GF DP -O ${1}.table
```
Load R
```
module load StdEnv/2020 java/13.0.2 r/4.1.2
```
Need to change file names twice in this script
```
library (ggplot2)
setwd("./")
# get all the files with the site data
files <- list.files(path = ".", pattern = "allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz.table")

my_data<-c()
temp<-c()

# read in the data and name the df based on the file name
for(f in 1:length(files)) {
  temp <- read.table(files[f], header = T)
  my_data <- rbind(my_data,temp)
  temp<-c()
}  

# my_data now has coverage per site info for each sample for all sites, genomewide
dim(my_data)
colnamez <- colnames(my_data)

# here is a function to calculate moving averages
# https://stackoverflow.com/questions/743812/calculating-moving-average
moving_fun <- function(x, w, FUN, ...) {
  # x: a double vector
  # w: the length of the window, i.e., the section of the vector selected to apply FUN
  # FUN: a function that takes a vector and return a summarize value, e.g., mean, sum, etc.
  # Given a double type vector apply a FUN over a moving window from left to the right, 
  #    when a window boundary is not a legal section, i.e. lower_bound and i (upper bound) 
  #    are not contained in the length of the vector, return a NA_real_
  if (w < 1) {
    stop("The length of the window 'w' must be greater than 0")
  }
  output <- x
  for (i in 1:length(x)) {
    # plus 1 because the index is inclusive with the upper_bound 'i'
    lower_bound <- i - w + 1
    if (lower_bound < 1) {
      output[i] <- NA_real_
    } else {
      output[i] <- FUN(x[lower_bound:i, ...])
    }
  }
  output
}

# example
# v <- seq(1:10)

# compute a MA(2)
# moving_fun(v, 2, mean)


# make a new dataframe that has the moving average for each sample throughout the 
# genome.  No worries if the window goes across chromosomes
mv_ave_df2 <- c()
for (i in 3:length(my_data)){  
  # calculate the moving average for each column
  moving_average <- moving_fun(my_data[,i], 50, mean)
  # add this to a new dataframe
  mv_ave_df2 <- cbind(mv_ave_df2,moving_average)
  # rename the column to match the sample
  colnames(mv_ave_df2)[i-2] = colnamez[i]
}  

# add chromosome and position data to mv_ave_df2
mv_ave_df3 <- data.frame(mv_ave_df2,my_data$CHROM,my_data$POS)

colnames(mv_ave_df3)[7] <- "CHROM"
colnames(mv_ave_df3)[8] <- "POS"


# calculate mean depth and sd per sample
mean_depth<- c()
sd_depth<- c()

for (i in 3:length(my_data)){
  x<- mean(my_data[,i],na.rm=T)
  mean_depth<-append(mean_depth,x, after = length(mean_depth))
  y<-sd(my_data[,i],na.rm=T)
  sd_depth<-append(sd_depth,y, after = length(sd_depth))
}  

cutoff_vector <- c()


# identify chr and positions in any sample that are >3 sd above that samples mean coverage
# based on the rolling average
# first make a vector with cutoff values for each sample
for (i in 3:ncol(mv_ave_df3)){
  cutoff <- mean_depth[i-2] + 3*sd_depth[i-2]
  cutoff_vector <- append(cutoff_vector,cutoff, after=length(cutoff_vector))
}

mean_over_cuttoff <- mean_depth/cutoff_vector
par(mfrow=c(2,2), mar=c(2,2,2,2)) 
plot(mean_over_cuttoff)
# give the elements in cutoff_vector some names
names(cutoff_vector) <- c('BJE4441_female_west_merge_sorted_dedup.bam.DP',
                          'BJE4442_male_west_merge_sorted_dedup.bam.DP',
                          'BJE4515_female_east_merge_sorted_dedup.bam.DP',
                          'BJE4536_male_east_merge_sorted_dedup.bam.DP',
                          'SRR6357672.1_trim_fixed_merged_dedup.bam.DP',
                          'SRR6357673_trim_sorted_dedup.bam.DP' )

# now cycle through each column and identify chr and pos of the bad ones
sub <- NULL
sub <- subset(mv_ave_df3, 
              BJE4441_female_west_merge_sorted_dedup.bam.DP > cutoff_vector["BJE4441_female_west_merge_sorted_dedup.bam.DP"]  | 
                BJE4442_male_west_merge_sorted_dedup.bam.DP > cutoff_vector["BJE4442_male_west_merge_sorted_dedup.bam.DP"] |
                BJE4515_female_east_merge_sorted_dedup.bam.DP  > cutoff_vector["BJE4515_female_east_merge_sorted_dedup.bam.DP"] |
                BJE4536_male_east_merge_sorted_dedup.bam.DP  > cutoff_vector["BJE4536_male_east_merge_sorted_dedup.bam.DP"] |
                SRR6357672.1_trim_fixed_merged_dedup.bam.DP  > cutoff_vector["SRR6357672.1_trim_fixed_merged_dedup.bam.DP"] |
                SRR6357673_trim_sorted_dedup.bam.DP > cutoff_vector["hecki_PF647.DP"] |
                SRR6357673_trim_sorted_dedup.bam.DP> cutoff_vector["SRR6357673_trim_sorted_dedup.bam.DP"] 
              )

dim(sub)
dim(mv_ave_df3)


to_file <- data.frame(sub$CHROM, sub$POS)

# write to file
write.table(to_file, "Chr8L_positions_to_exclude.txt", append = FALSE, sep = "\t", dec = ".",
            row.names = F, col.names = F)
```

Remove quotes
```
sed -i 's/"//g' Chr8L_positions_to_exclude.txt
```
Now remove these sites
```
module load nixpkgs/16.09 intel/2018.3 vcftools/0.1.16
vcftools --gzvcf allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz --out allsites_Chr8L_SNPs.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz.recode.vcf.gz --exclude-positions Chr8L_positions_to_exclude.txt --recode
```
