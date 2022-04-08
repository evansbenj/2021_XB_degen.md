# Lumpy analysis

Following suggested workflow here:
https://github.com/arq5x/lumpy-sv

# Align specifically for lumpy analysis

```
#!/bin/sh
#SBATCH --job-name=bwa_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=120:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa_align.%J.out
#SBATCH --error=bwa_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2020_align_paired_fq_to_ref.sh pathandname_of_ref path_to_paired_fq_filez
# sbatch 2020_align_paired_fq_to_ref.sh /home/ben/projects/rrg-ben/ben/2018_Austin_XB_genome/Austin_genome/Xbo.v1.fa.gz
 pathtofqfilez

module load bwa/0.7.17
module load samtools/1.10
module load StdEnv/2020 samblaster/0.1.26

for file in ${2}/*_trim.R1.fq.gz ; do         # Use ./* ... NEVER bare *    
    if [ -e "$file" ] ; then   # Check whether file exists.
	bwa mem ${1} ${file::-14}_trim.R1.fq.gz ${file::-14}_trim.R2.fq.gz -t 16 |  samblaster --excludeDups --addMateT
ags --maxSplitCount 2 --minNonOverlap 20 | samtools view -S -b - > ${file::-14}_lumpy.bam
    fi
done
```
# Merge if multiple fastq files aligned separately
```
#!/bin/sh
#SBATCH --job-name=merge
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=32:00:00
#SBATCH --mem=16gb
#SBATCH --output=merge.%J.out
#SBATCH --error=merge.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2021_samtools_merge4.sh merged_path_and_file in1_path_and_file in2_path_and_file in3_path_and_file in4_path_an
d_file
module load samtools/1.10

samtools merge ${1} ${2} ${3} ${4} ${5}
samtools index ${1)
```
# Extract discordants and splitters
```
#!/bin/sh
#SBATCH --job-name=discordant
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=32:00:00
#SBATCH --mem=16gb
#SBATCH --output=discordant.%J.out
#SBATCH --error=discordant.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2022_samtools_extractdiscordant.sh path_and_file_bam
module load samtools/1.10

# Extract the discordant paired-end alignments.
 samtools view -b -F 1294 ${1} > ${1}.discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h ${1} \
    | /home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/bin/lumpy-sv/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > ${1}.splitters.unsorted.bam

# Sort both alignments
samtools sort ${1}.discordants.unsorted.bam -o ${1}.discordants.sorted.bam
samtools sort ${1}.splitters.unsorted.bam -o ${1}.splitters.sorted.bam
```
# Run lumpy on multiple files
```
#!/bin/sh
#SBATCH --job-name=lumpyexpress
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=32:00:00
#SBATCH --mem=16gb
#SBATCH --output=lumpyexpress.%J.out
#SBATCH --error=lumpyexpress.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2022_samtools_extractdiscordant.sh path_and_file_bam
module load python/3.8
module load StdEnv/2020
module load scipy-stack/2020b
module load samtools/1.10
module load samblaster/0.1.26 sambamba/0.8.0 lumpy/0.2.13

lumpyexpress \
    -K /home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/ben_scripts/lumpyexpress.config \
    -B sample1.bam,sample2.bam,sample3.bam \
    -S sample1.splitters.bam,sample2.splitters.bam,sample3.splitters.bam \
    -D sample1.discordants.bam,sample2.discordants.bam,sample3.discordants.bam \
    -o multi_sample.vcf
```

# Smoove

Smoove is a new wrapper for lumpy that implements filtering steps and streamlines genotyping of structural variants:
https://brentp.github.io/post/smoove/

The files need to be named in a specific way that matches their readgroups.  I did this here on graham (all are sorted and indexed):
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/lumpy_linkz
```
There are lots of dependencies for smoove.  Some are available as modules on computecanada and I installed the others:
```
module load StdEnv/2020 python/2.7.18
pip install git+https://github.com/hall-lab/svtyper.git
had to do a download and local install of a python module to get svtools to install:
python -m pip install ./pandas-0.19.2.tar.gz 
```
had to tweak something for git to get mosdepth to install:
```
git config --global url."https://github.com/".insteadOf git://github.com/
```
I added lots of paths to my bash_profile to make them run smoothly
```
emacs ~/.bash_profile
smoove call -x --name my-cohort --exclude $bed --fasta $reference_fasta -p $threads --genotype /path/to/*.bam
smoove call -x --name testipooh /home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/Austin_genome/Xbo.v1.fa --genotype /home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/lumpy_linkz/*.bam
```

# Visualization with Samplot
Useful video to demonstrate capability:
```
https://www.youtube.com/watch?v=-2Iyq3gkXGM
```
Installation instead of conda:
```
module load StdEnv/2020 python/3.10.2
pip install samplot
pip install pysam
```

# Samplot
I got this to work locally:
```
samplot plot -n BJE4515 -b WGS_chr8L_7S_bams/BJE4515_female.split_Chr7S.bam -o BJE4515_Chr7S_10_21Mb_split.png -c Chr7S -s 10000000 -e 21000000 -t INV --maxmb 20
```
But not yet on computecanada.


# BELOW NOT USED

# Analyze bam with lumpy
```
#!/bin/sh
#SBATCH --job-name=lumpy
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=24gb
#SBATCH --output=lumpy.%J.out
#SBATCH --error=lumpy.%J.err
#SBATCH --account=def-ben

# run like this:
# sbatch Lumpy.sh bam_path_and_bam_prefix_without_thebamsuffix

lumpy \
    -mw 4 \
    -tt 0 \
    -pe id:BJE4441,bam_file:BJE4441_female_west_merge_sorted_dedup.bam_rg.bam_Chr7S.discordants.bam,histo
_file:BJE4441_female_west_merge_sorted_dedup.bam_rg.bam_Chr7S.lib1.histo,mean:303.7,stdev:118.4,read_leng
th:151,min_non_overlap:151,discordant_z:5,back_distance:10,weight:1,min_mapping_threshold:20 \
    -sr id:BJE4441,bam_file:BJE4441_female_west_merge_sorted_dedup.bam_rg.bam_Chr7S.splitters.bam,back_di
stance:10,weight:1,min_mapping_threshold:20 \
    > BJE4441_female_west_merge_sorted_dedup.bam_rg.bam_Chr7S.vcf
```

# sort the vcf file by coordinate
```
cat BJE4441_female_west_merge_sorted_dedup.bam_rg.bam_Chr7S.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' > BJE4441_female_west_merge_sorted_dedup.bam_rg.bam_Chr7S_sorted.vcf
```
# compress and index using tabix
```
bgzip -c BJE4441_female_west_merge_sorted_dedup.bam_rg.bam_Chr7S_sorted.vcf > BJE4441_female_west_merge_sorted_dedup.bam_rg.bam_Chr7S_sorted.vcf.gz
tabix -p vcf BJE4441_female_west_merge_sorted_dedup.bam_rg.bam_Chr7S_sorted.vcf.gz
```
# Use lumpy to identify structural variants

# Making a json file for each bam file (this can also be done with svtyper I think - check readme)

In this directory:
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/XB_wild_east_west_WGS_data/bam_dedup_rg
```
```
module load gcc bwa samtools bcftools blast muscle primer3 python/2.7
virtualenv --no-download vcf_env
source vcf_env/bin/activate
pip install git+https://github.com/hall-lab/svtyper.git
pip install gdc-readgroups --user
pip install cwltool --user
wget https://raw.githubusercontent.com/NCI-GDC/gdc-readgroups/master/Dockerfile
wget https://raw.githubusercontent.com/NCI-GDC/gdc-readgroups/master/gdc-readgroups.cwl
cwltool gdc-readgroups.cwl --INPUT BJE4536_male_east_merge_sorted_dedup.bam_rg.bam
gdc-readgroups bam-mode --b BJE4536_male_east_merge_sorted_dedup.bam_rg_Chr7S.bam
gdc-readgroups bam-mode --b BJE4515_female_east_merge_sorted_dedup.bam_rg_Chr7S.bam
gdc-readgroups bam-mode --b BJE4442_male_west_merge_sorted_dedup.bam_rg.bam_Chr7S.bam
gdc-readgroups bam-mode --b BJE4441_female_west_merge_sorted_dedup.bam_rg.bam
```
# Analyze the lumpy vcf file with sctyper
Directory:
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/XB_wild_east_west_WGS_data/bam_dedup_rg/svtyper
```
```
./svtyper.sh ../BJE4515_female_east_merge_sorted_dedup.bam_rg_Chr7S_sorted ../BJE4515_female_east_merge_sorted_dedup.bam_rg_Chr7S
```

# Analyze bam with ClipSV
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/XB_wild_east_west_WGS_data/bam_dedup_rg/ClipSV
```

```
module load StdEnv/2020 minimap2/2.24
module load samtools bedtools
module load velvet
source activate python3
PATH=$PWD/ClipSV/:$PATH
python3 ./clipsv.py -t 1 -b ../BJE4515_female_east_merge_sorted_dedup.bam_rg_Chr7S.bam -g ../../../Austin_genome/Xbo_Chr7S.fa
```

