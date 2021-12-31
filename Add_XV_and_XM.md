# Adding XV and XM sample for Chr7S

path to XM sample:
```
/home/ben/projects/rrg-ben/ben/2021_XLXGXM/raw_data_not_demultiplexed/samples/concat/XM_1.fq.gz
```
path to XV sample
```

```

Align
```
#!/bin/sh
#SBATCH --job-name=bwa_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=48:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa_align.%J.out
#SBATCH --error=bwa_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2020_align_paired_fq_to_ref.sh pathandname_of_ref path_and_name_singleend_fq_filez


module load bwa/0.7.17
module load samtools/1.10

echo ${file::-6}
echo bwa mem ${1} ${file::-6}.fq.gz -t 16 | samtools view -Shu - | samtools sort - -o ${file::-6}_sorted.bam
bwa mem ${1} ${file::-6}.fq.gz -t 16 | samtools view -Shu - | samtools sort - -o ${file::-6}_sorted.bam
samtools index ${file::-6}_sorted.bam
```

