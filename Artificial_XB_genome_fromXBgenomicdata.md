# Making fake XB genome from mapped XB reads to XL genome

Working in this directory:
```
/home/ben/projects/rrg-ben/ben/2021_Austin_XB_genome/MandF_WGS_filtered_and_removed_vcfs_XL_ref
```
Use bcftools to make a consensus file from mapped reads:
```
#!/bin/sh
#SBATCH --job-name=bcftools
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --mem=2gb
#SBATCH --output=bcftools.%J.out
#SBATCH --error=bcftools.%J.err
#SBATCH --account=def-ben

# This program will make a consensus out of a vcf file from XB using the genome of XL.
# It will insert IUPAC symbols for any variant and replace XL positions with homoz XB postions too
# It will do this using variants from the female genome (SRR6357673) that are mapped to the XL genome

# execute like this: ./2021_bcftools_extract_sections_from_vcf.sh path_and_filename_of_vcf.gz chr
# load these modules before running:
module load StdEnv/2020 gcc/9.3.0 bcftools/1.11
module load StdEnv/2020  gcc/9.3.0 samtools/1.13

samtools faidx /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_refgenome/XENLA_9.2_genome.fa ${2} | bcftools consensus -H A -s SRR635
7673_trim._sorted.bam ${1} > ${1}_XLmapping_consensus.fa
```
I had to use the '-H A' flag because I had issues with the '-I' flag for chr1L only for some reason.

After concatenating the files, I had to replace asterisks with an 'N':
```
sed -i 's/\*/N/g' ../Artifical_XB_genome/artificial_reference/SRR6357673_artificial_XBgenome_basedonmappingtoXL.fa
```
because haplotype caller didn't like ascii char 42 in the ref. 
