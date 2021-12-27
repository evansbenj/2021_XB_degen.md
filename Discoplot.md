# Plotting structural variation with Discoplot

I used git locally to install Discoplot from git. On graham, I then installed dependencies like this:
```
pip instal -r requirements.txt
```
I changed to executable:
```

```
And tried this command:
```
./DiscoPlot.py -bam BJE4441_female_west_merge_sorted_dedup.bam_rg.bam -bin 10000 -o discoplot.png
```
Which is running now.

I'm also subsetting the bam file to include only chr7 like this:
```
samtools view -b BJE4515_female_east_merge_sorted_dedup.bam_rg.bam Chr7S > BJE4515_female_east_merge_sorted_dedup.bam_rg_Chr7S.bam
```
