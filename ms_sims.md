# `ms` complilation
Didn't work at first but this worked for a local install:
```
arch -arch x86_64 gcc -O3 -o ms ms.c streec.c rc -lmc
```
# `ms` commandline
```
./ms 6 100 -t 10 -I 3 2 2 2 -n 1 0.125 -n 2 0.375 -n 3 0.500 -ej 0.1 1 3 -ej 0.1 2 3
```
This generates 100 simulations with:
* sample size = 6
* theta = 10
* `-I 3 2 2 2` : 2 samples from each of 3 populations; no migration
* `-n 1 0.125 -n 2 0.375 -n 2 0.500`: pop1 to be 1/8th of N_ancestral (W-chr), pop2 to be 3/8ths of N_ancestral (Z-chr) and pop3 to be 50% of N_ancestral (west pop)
* `-ej 0.1 2 3 -ej 0.1 2 3`: pop1 and pop2 coalesce with pop3 at time 0.1

From this I can sample that are heteroz in artificial genotypes from 1 chr from pop1 (Wchr) and 1 chr from pop2 (Zchr) (2 females) and homozy in artifical genotypes from 1 chr from pop1 (Wchr) and 1 chr from pop2 (Zchr) (2 males) and then calculate genotype frequencies of these sites in pop3. 

This could be done also with a range of divergence times (the first parameter currently set to 01 in the `-ej` command)

# Playing with one simulation:
```
./ms 6 1 -t 100 -I 3 2 2 2 -n 1 0.125 -n 2 0.375 -n 3 0.500 -ej 0.1 1 3 -ej 0.1 2 3
```

I increased theta to be 100 to generate lots of polymorphism on 1 chr.  This could be increased to a large number - each site has an independent coalesent history anyhow so this may be similar to doing independent simulations and sampling only 1 variant.  Need to look into this.
