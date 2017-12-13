# bottleneck
A forward genetic simulator for playing around with demographic bottlenecks

This repository contains C code for a forward-time genetic simulation program, two simple Perl scripts that provide sample values for two scenarios, a data file containing the derived allele frequency distribution for the African superpopulation in the 1000 Genomes Project, and some R code that can be used to display the results of simulations.

The demographic model consists of 3 eras:
- An initial constant-sized population (simulated by a burn-in period of 10 x (initial population size) generations). 
- A first era of exponential growth.
- A second era of exponential growth.
A bottleneck in population size can be imposed before either the first or second growth eras, and is specified by the new, shrunken population size. Subsequent growth will start from this size.

The population is modeled as 2 x pop-size independent chromosomes. Population growth is deterministic and discrete -- which means that specifying a growth rate of 0.1% for a population of size 100 will not do anything, since the population cannot grow by 0.2 chromosomes per generation. 

One can also specify a selection coefficient to apply purifying selection (the original purpose of this code).

## Installation

The C source code must be compiled before use, along with associated functions (supplied with package). Using the current Mac  built-in C compiler, I compile it with the command

```
cc -O3 -lm -o forward -Wall forward.c brano_x.c prano_x.c gammln.c random_x.c 
```
which should yield the executable file *forward*, but other compilers should work as well.

## Execution

*forward* is run from the command line. Model parameters are supplied through a bevy of command line arguments when executing the program. The user must specify, in order, the following:


- selection coefficient (positive for purifying selection); set to 0.0 for neutral drift.
- initial population size (number of diploid individuals).
- fractional growth rate for first growth era; a value of 0.01 corresponds to 1% growth per generation.
- length of first growth era in generations.
- bottleneck population size.
- fractional growth rate for 2nd growth era.
- length of second growth era in generations.
- optional sample size, if final output represents a small sample of full population. (I haven't tested this option in a long time.)
- mutation rate, in new mutations per generation. One generally simulates a tiny genome (e.g. 0.5 mutations/generation, or about 1% of the true rate for humans) and scales the output up.
- starting time for purifying selection to kick in, in generations (after burn-in); set to -1 for constant selection. Irrelevant if s = 0.
- output file name; output files will have ".txt" and ".log" appended.
- flag to indicate when bottleneck should occur. 0 means before first growth era, 1 means before second era.

The allele frequency distribution output will be in file <filename>.txt and a log of the run in <filename>.log.

The provided sample scripts, run_250kya_8k.pl and run_100kya_4k.pl, show the parameters I used to generate scenarios with a bottleneck 250,000 years ago/final population size ~8000 and with a bottleneck 100,000 years ago/final population size ~4000, respectively. 

I have included some R code. I use it to scale up simulated output to compare with real human genomes, normalize the pre-bottleneck population size to match the 1000 Geneomes data in the 0.6 to 0.7 range, and plot both simulated and real data. If you set the working directory in R to be where you have put the code and 
```
source("plotPublic.R")
```
it should plot the 100ky/4k data (assuming you have the appropriate output file lying about).
