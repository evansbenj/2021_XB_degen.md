# Plotting structural variation with Discoplot

I used git locally to install Discoplot from git. On graham, I then installed dependencies like this:
```
pip install -r requirements.txt
```
The above did not work. But this maybe will:
```
module load python/2.7
virtualenv /home/$USER/my_venv
source /home/$USER/my_venv/bin/activate
pip install numpy
pip install matplotlib
pip install pysam
./DiscoPlot.py -bam ../DiscoPlot.py -bam ../BJE4536_male_east_merge_sorted_dedup.bam_rg_Chr7S.bam -bin 10000 -dpi 200 -o temp_Chr7S_10000_disco.png
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
