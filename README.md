# bottleneck
A forward genetic simulator for playing around with demographic bottlenecks

This repository contains C code for a forward-time genetic simulation program, two simple Perl scripts that provide sample values for two scenarios, a data file containing the derived allele frequency distribution for the African superpopulation in the 1000 Genomes Project, and some R code that can be used to display the results of simulations.

## Installation

The C source code must be compiled before use, along with associated functions (supplied with package). Using the current Mac  built-in C compiler, I compile it with the command

```
cc -O3 -lm -o forward -Wall forward.c brano_x.c prano_x.c gammln.c random_x.c 
```

which should yield the executable file *forward*, but other compilers should work as well.

## Execution

*forward* is run from the command line. Model parameters are supplied through a bevy of command line arguments when executing the program


- -i: File of genotype data. See below for format.
