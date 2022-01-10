# `ms` complilation
Didn't work at first but this worked for a local install:
```
arch -arch x86_64 gcc -O3 -o ms ms.c streec.c rc -lmc
```
# `ms` commandando for Scenario 2
```
#!/usr/bin/env perl
use strict;
use warnings;

# This program does a bunch of simulations for Xborealis sex chromosomes
# using ms and then summarizes the results

my $population_divergence_limit=0.5;
my $x;
my $y;
my $status;

# We need to simulate over a range of divergence values for tao0 (pop divergence)
# and tao1 (recomb suppression)
for ($x = 0 ; $x <= $population_divergence_limit ; $x+=0.01) {
	for ($y = 0 ; $y <= $x ; $y+=0.01) {
		print "Tao0 ",$x," tao1 ",$y,"\n";
		# perform the first simulation
		$status = system("./msdir/ms 12 100000 -t 100 -I 3 2 6 4 -n 1 0.083 -n 2 0.250 -n 3 0.667 -ej ".$y." 1 2 -ej ".$x." 2 3 > out");
		# now summarize the simulations
		$status = system("./XB_SNPcount_EastWest_ms_sims_new.pl out out2");
	}	# end of $y
} # end of $x			
```
This generates 100 simulations with:
* sample size = 12 chromosome (6 diploids)
* theta = 4Neu = 100
* `-I 3 2 6 4` : 2 samples of the W chr from pop1, 6 samples of the Z chr from pop2, 4 samples from pop3; no migration
* `-n 1 0.083 -n 2 0.250 -n 3 0.667`: pop1 to be 1/12th of N_ancestral (W-chr), pop2 to be 3/12ths of N_ancestral (Z-chr) and pop3 to be 2/3rds of N_ancestral (west pop)
* `-ej ".$y." 1 2 -ej ".$x." 2 3`: pop1 and pop2 coalesce at time $y and then these coalesce with pop3 at time $x

# `ms` commandando for Scenario 3a
```
#!/usr/bin/env perl
use strict;
use warnings;

# This program does a bunch of simulations for Xborealis sex chromosomes
# using ms and then summarizes the results

my $population_divergence_limit=0.5;
my $x;
my $y;
my $status;

# We need to simulate over a range of divergence values for tao0 (pop divergence)
# and tao1 (recomb suppression)
for ($x = 0 ; $x <= $population_divergence_limit ; $x+=0.01) {
	for ($y = 0 ; $y <= $x ; $y+=0.01) {
		print "Tao0 ",$x," tao1 ",$y,"\n";
		# perform the first simulation
		# pop1 is the east W; pop2 is the eastZ; pop3 is the west
		# Scenario3a: Z becomes an autosome in west
		# Scenario3b: W becomes an autosome in west
		# this is Scenario3a:
		$status = system("./msdir/ms 12 100000 -t 100 -I 3 2 6 4 -n 1 0.083 -n 2 0.250 -n 3 0.667 -ej ".$y." 3 2 -ej ".$x." 1 2 -en ".$y." 2 0.273 > out3a");
		# now summarize the simulations
		$status = system("./XB_SNPcount_EastWest_ms_sims_new.pl out3a out_all_3a");
	}	# end of $y
} # end of $x			
```
# `ms` commandando for Scenario 3b
```
#!/usr/bin/env perl
use strict;
use warnings;

# This program does a bunch of simulations for Xborealis sex chromosomes
# using ms and then summarizes the results

my $population_divergence_limit=0.5;
my $x;
my $y;
my $status;

# We need to simulate over a range of divergence values for tao0 (pop divergence)
# and tao1 (recomb suppression)
for ($x = 0 ; $x <= $population_divergence_limit ; $x+=0.01) {
	for ($y = 0 ; $y <= $x ; $y+=0.01) {
		print "Tao0 ",$x," tao1 ",$y,"\n";
		# perform the first simulation
		# pop1 is the east W; pop2 is the eastZ; pop3 is the west
		# Scenario3a: Z becomes an autosome in west
		# Scenario3b: W becomes an autosome in west
		# this is Scenario3b:
		$status = system("./msdir/ms 12 100000 -t 100 -I 3 2 6 4 -n 1 0.083 -n 2 0.250 -n 3 0.667 -ej ".$y." 3 1 -ej ".$x." 1 2 -en ".$y." 1 0.111 > out3b");
		# now summarize the simulations
		$status = system("./XB_SNPcount_EastWest_ms_sims_new.pl out3b out_all_3b");
	}	# end of $y
} # end of $x			
```

# Summarizing sims
```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use List::Util qw( min max );


#  This program reads in the output from ms sims.

# For each simulation, it will first identify eligible sites, and from these
# it will count the number of double homoz W-specific
# double homoz Z-specific, and double hets in the west.

# it will then quantify the proportion of simulations where these two ratios:
# double homoz Z-specific:double homoz W-specific
#        and
# double homoz Z-specific:double hets

# are both within 10% of the observed.

# the first two samples are the east W.  The next 6 are the east Z
# the last 4 are the two diploid from the west.

# here is an example commandline:
# ././msdir/ms 12 1000000 -t 100 -I 3 2 6 4 -n 1 0.083 -n 2 0.250 -n 3 0.667 -ej tao1 1 2 -ej tao0 2 3 > out
# where tao0 is the time of population divergence and tao1 is the time of recombination suppression in the east

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];


unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file 1.\n";
	exit;
}


unless (open(OUTFILE, ">>$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
#print "Creating output file: $outputfile\n";

my $pattern="";
my @temp;

# in order to be accepted, the pattern must begin either with this:
# 11000000
# or this:
# 00111111
# because this means the two east females are hets and the two east males are homoz
# the first east female is allele 1+3, the second east female is allele 2+4
# the first east male is allele 5+6 and the second east male is allele 7+8 
# the last simulation will be ignored

my $double_homoz_Wspecific=0;
my $double_homoz_Zspecific=0;
my $double_hets=0;
my $firstline=0;
my %block=();
my $n;
my $m;
my $individ=0;
my $individuals=12;
my $sites=0;
my $acceptable=0;
my $accepted=0;
my $observed_ratio_1 = 4.9;
my $observed_ratio_2 = 38.3;
my $percentage = 0.1;
my $sims=0;
my $sum_accepted_proportions=0;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	if($firstline==2){
		@temp=split(//,$line);
		if(defined($temp[0])){
			if(($temp[0] eq 0)||($temp[0] eq 1)){					
				for ($n = 0 ; $n <= $#temp ; $n+=1){
					$block{$n}[$individ]=$temp[$n];
				}
				$individ+=1;
				$sites=$#temp;
			}
			# now we have a block hash; each site pattern is $block{$n}[$individ=0-11]
			elsif(($temp[0] eq '/')&&($sites > 0)){ # check if the pattern is ok, and if so characterize it
				# cycle through all sites
				for ($n = 0; $n <= $sites; $n+=1){
					$pattern="";
					# make a pattern
					for ($m = 0; $m < $individuals; $m+=1){
						$pattern = $pattern.$block{$n}[$m];
					}				
					if(($pattern =~ /^11000000/)||($pattern =~ /^00111111/)){
						$acceptable+=1;
						# count double homoz for W-specific
						if(($pattern eq '110000001111')||($pattern eq '001111110000')){
							$double_homoz_Wspecific+=1;
							#print "double_homoz_Wspecific ",$pattern,"\n";
						}
						elsif(($pattern eq '110000000000')||($pattern eq '001111111111')){
							$double_homoz_Zspecific+=1;
							#print "double_homoz_Zspecific ",$pattern,"\n";
						}
						elsif(($pattern eq '110000001010')||($pattern eq '110000000101')||
						($pattern eq '110000001001')||($pattern eq '110000000110')||
						($pattern eq '001111111010')||($pattern eq '001111110101')||
						($pattern eq '001111111001')||($pattern eq '001111110110')){
							$double_hets+=1;
						}
					}
				}
				if(($double_homoz_Wspecific != 0)&&($double_hets != 0)){
					if(
					($double_homoz_Zspecific/$double_homoz_Wspecific > ($observed_ratio_1-$percentage*$observed_ratio_1))&&
					($double_homoz_Zspecific/$double_homoz_Wspecific < ($observed_ratio_1+$percentage*$observed_ratio_1))&&
					($double_homoz_Zspecific/$double_hets > ($observed_ratio_2-$percentage*$observed_ratio_2))&&
					($double_homoz_Zspecific/$double_hets < ($observed_ratio_2+$percentage*$observed_ratio_2))){
						$accepted+=1;
					}
				}	
				%block=();
				$individ=0;
				$sites=0;
				$double_homoz_Wspecific=0;	
				$double_homoz_Zspecific=0;
				$double_hets=0;
				# weight the proportion equally for each simulation
				$sims+=1;
				if($acceptable>0){
					print $accepted/$acceptable,"\n";
					$sum_accepted_proportions+=$accepted/$acceptable; # this will be divided by the number of sims below
				}
				# reset counters for each sim	
				$accepted=0;
				$acceptable=0;
			}
		}
	}
	elsif($firstline==1){
		$firstline=2;
	}
	elsif($firstline==0){
		print OUTFILE $line,"\t";
		$firstline=1;
	}
	
} # end while

print "accepted ratio: ",$sum_accepted_proportions/$sims," eligible_sites ",$acceptable,"\n";
print OUTFILE "mean_accepted_ratio: ",$sum_accepted_proportions/$sims," eligible_sites ",$acceptable,"\n";
```
## Contour plot
```R
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2021_Xborealis_sexchr_degen/ms")
par(mfrow=c(2,2), mar=c(2,2,2,2)) 
simulations <- read.table("ms_out.txt", header = T)
t0 <- simulations$T0[simulations$T0<0.40]
t1 <- simulations$T1[simulations$T0<0.40]
lnL <- simulations$lnL[simulations$T0<0.40]
#View(simulations)

lnL_matrix <- matrix(nrow=40,ncol=40)
simulations[832,"T0"]
# make a new matrix out of the dat df
for(i in seq(from=1, to=nrow(simulations))){
  #print(paste(100*simulations[i,"T0"]+1,100*simulations[i,"T1"]+1,-simulations[i,"lnL"]))
  lnL_matrix[100*simulations[i,"T0"]+1,100*simulations[i,"T1"]+1]<- -simulations[i,"lnL"]
}
lnL_matrix[29,1]<- 0
View(lnL_matrix)
x<-seq(0,0.39,0.01) # this will be t0
y<-seq(0,0.39,0.01) # this will be t1
#jet.colors <-colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
png(filename = "contour.png",w=500, h=500,units = "px", bg="transparent")
  filled.contour(x, y, lnL_matrix,
                 nlevels = 100,
               lwd = 2, lty = 1,
               color.palette=colorRampPalette(c('white','white','white','white','white','white','white','white','white','white','white','light blue','blue','yellow','red','darkred','black')),
               #color.palette=colorRampPalette(c('white','blue','yellow','red','darkred')),
               #olor = jet.colors, 
               #color.palette = colorRampPalette(YlOrBr, space = "Lab"),
               #color.palette = colorRampPalette(c("red", "white", "blue")),
               #color = terrain.colors,
               xlab = "t0", 
               ylab = "t1",
               key.title = {par(cex.main=1);title(main="        -lnL")},
               )
dev.off()

```
## Simulating panmixia (applies to the recombining portion of Chr8L; not used)
```
./msdir/ms 12 100000 -t 1 -I 2 8 4 -n 1 0.5 -n 2 0.500 -ej .2 1 2 > out
```
Here is a perl script to parse the ms sims:
```
#!/usr/bin/env perl
use strict;
use warnings;
use List::MoreUtils qw/ uniq /;
use List::Util qw( min max );


#  This program reads in the output from ms sims.

# It will take the first variable position from each simulation,
# construct genotypes, and check the number of double homoz W-specific
# double homoz Z-specific, and double hets in the east.

# the first two samples are the east W.  The next 6 are the east Z
# the last 4 are the two diploid from the west.

# here is an example commandline:
# ./ms 12 10 -t 10 -I 2 2 6 4 -n 1 0.125 -n 2 0.375 -n 3 0.500 -ej .42 1 3 -ej .42 2 3 > out 

my $inputfile = $ARGV[0];
my $outputfile = $ARGV[1];


unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file 1.\n";
	exit;
}


unless (open(OUTFILE, ">>$outputfile"))  {
	print "I can\'t write to $outputfile\n";
	exit;
}
#print "Creating output file: $outputfile\n";

my $pattern="";
my @temp;

# in order to be accepted, the pattern must begin either with this:
# 11000000
# or this:
# 00111111
# because this means the two east females are hets and the two east males are homoz
# the first east female is allele 1+3, the second east female is allele 2+4
# the first east male is allele 5+6 and the second east male is allele 7+8 
# the last simulation will be ignored

my $double_homoz_Wspecific=0;
my $double_homoz_Zspecific=0;
my $double_hets=0;
my $firstline=0;

while ( my $line = <DATAINPUT>) {
	chomp($line);
	if($firstline==0){
		print OUTFILE $line,"\n";
		$firstline=1;
	}
	@temp=split(//,$line);
	if(defined($temp[0])){
		if(($temp[0] eq 0)||($temp[0] eq 1)){
			$pattern = $pattern.$temp[0];
		}
		elsif($temp[0] eq '/'){ # check if the pattern is ok, and if so characterize it
			if(($pattern =~ /^11000000/)||($pattern =~ /^00111111/)){
				#print $pattern,"\n";
				# count double homoz for W-specific
				if(($pattern eq '110000001111')||($pattern eq '001111110000')){
					$double_homoz_Wspecific+=1;
				}
				elsif(($pattern eq '110000000000')||($pattern eq '001111111111')){
					$double_homoz_Zspecific+=1;
				}
				elsif(($pattern eq '110000001010')||($pattern eq '110000000101')||
				($pattern eq '110000001001')||($pattern eq '110000000110')||
				($pattern eq '001111111010')||($pattern eq '001111110101')||
				($pattern eq '001111111001')||($pattern eq '001111110110')){
					$double_hets+=1;
				}
			}
			$pattern="";	
		}	
	}
} # end while
print "Zdoublehomoz/Wdoublehomoz: ",$double_homoz_Zspecific/$double_homoz_Wspecific,"\n";
print OUTFILE $double_homoz_Wspecific,"\t",$double_homoz_Zspecific,"\t",$double_hets,"\n";
```
