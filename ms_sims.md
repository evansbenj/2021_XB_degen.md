# `ms` complilation
Didn't work at first but this worked for a local install:
```
arch -arch x86_64 gcc -O3 -o ms ms.c streec.c rc -lmc
```
# `ms` commandline
```
./ms 12 100 -t 10 -I 3 2 6 4 -n 1 0.125 -n 2 0.375 -n 3 0.500 -ej 0.1 1 3 -ej 0.1 2 3
```
This generates 100 simulations with:
* sample size = 12 chromosome (6 diploids)
* theta = 10
* `-I 3 2 6 4` : 2 samples of the W chr from pop1, 6 samples of the Z chr from pop2, 4 samples from pop3; no migration
* `-n 1 0.125 -n 2 0.375 -n 2 0.500`: pop1 to be 1/8th of N_ancestral (W-chr), pop2 to be 3/8ths of N_ancestral (Z-chr) and pop3 to be 50% of N_ancestral (west pop)
* `-ej 0.1 2 3 -ej 0.1 2 3`: pop1 and pop2 coalesce with pop3 at time 0.1

From this I can sample that are heteroz in artificial genotypes from 1 chr from pop1 (Wchr) and 1 chr from pop2 (Zchr) (2 females) and homozy in artifical genotypes from 1 chr from pop1 (Wchr) and 1 chr from pop2 (Zchr) (2 males) and then calculate genotype frequencies of these sites in pop3. 

This could be done also with a range of divergence times (the first parameter currently set to 01 in the `-ej` command)

# Simulations
## Simulating the onset of recombination suppression from a panmictic ancestor in the east (applies to noon-recombining portion of Chr8L) 
This causes the population to split into 3 groups - the east W, the east Z and the west, with the relative proportions of Ne being 0.125, 0.375, and 0.5.
```
./msdir/ms 12 100000 -t 1 -I 3 2 6 4 -n 1 0.125 -n 2 0.375 -n 3 0.500 -ej .09 1 3 -ej .09 2 3 > out
```
But this still produces too many double hets.  I think this can be solved by introducing population structure where east and west diverge for a period of time before speciation of the sex chromosomes. I also am going to let the west population be 2X bigger than the east based on the observed polymorphism in non sex-linked regions

Here is the general structure of the sims:
./msdir/ms 12 1000000 -t 1 -I 3 2 6 4 -n 1 0.08325 -n 2 0.25005 -n 3 0.666 -ej tao1 1 2 -ej tao0 2 3 > out

Where tao1 is the time that recombination stops and tao0 is the time of divergence of east and west
and
tao1 <= tao0


need to modify stuff below

These two work pretty well
```
./msdir/ms 12 1000000 -t 1 -I 3 2 6 4 -n 1 0.125 -n 2 0.375 -n 3 0.500 -ej .05 1 2 -ej .11 2 3 > out
./msdir/ms 12 10-0000 -t 1 -I 3 2 6 4 -n 1 0.125 -n 2 0.375 -n 3 0.500 -ej .03 1 2 -ej .2 2 3 > out
```

## Simulating panmixia (applies to the recombining portion of Chr8L)
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
