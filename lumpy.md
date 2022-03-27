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
