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
