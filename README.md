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


- selection coefficient (positive float
