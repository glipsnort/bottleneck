#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "random.h"
#include "brano.h"
#include "prano.h"

void next_gen(double r_grow, double sel, int gen);

int popsize, nvar, *count, max_site, *birthgen;
double mut_rate;    // per single genome copy per generation
long unsigned int totmut=0;

int main(int argc, char **argv) {
  const int ndaf = 100;    // Number of bins to divide derived allele frequency distribution into for output
  double daf_binsize = 1. / ndaf;
  int nburn, startsize, expgen1, expgen2, sampsize, size, bottlesize, sel_start, dafbin;
  double sel, r_grow1, r_grow2, freq;
  double sel_store, r_grow, pi;
  int gen, ivar, iclass, sampled, nvarold, ntot_burn;
  int *count_by_class=NULL, *pre_by_class=NULL; 
  int *raw_distr=NULL, when_bneck, ntot_bneck=0, *raw_pre=NULL;
  FILE *outf=NULL, *logf=NULL;
  char filename[256], logfile[256];
  unsigned long seed;

  raw_distr = malloc(ndaf * sizeof(int));   // count of derived alleles binned by derived freq
  raw_pre = malloc(ndaf * sizeof(int));     // same, for variants predating bottleneck
  if (argc != 13) {
    fprintf(stderr, 
	    "Usage: forward <sel coeff (s>=0)> <initial pop size> <growth rate r1> <N gen expansion 1> \n");
    fprintf(stderr, "<bottleneck pop (neg=no bottlneck)> <growth rate r2> <N gen exp 2> \n");
    fprintf(stderr, "<sample size (neg=no sample), n diploids> <mutation rate> <start selection (gen)> ");
    fprintf(stderr, "<output file base> <bneck when [0=after burnin, 1=after 1st expansion]>\n");
    abort();
  }
  strcpy(logfile, argv[11]);
  strcat(logfile, ".log");
  logf = fopen(logfile, "w");
  if (logf == NULL) {fprintf(logf, "Could not open file %s\n", logfile); abort();}
  sel_store = strtod(argv[1], NULL);
  startsize = strtol(argv[2], NULL, 10);
  r_grow1 = strtod(argv[3], NULL);
  expgen1 = strtol(argv[4], NULL, 10);
  bottlesize = strtol(argv[5], NULL, 10);
  r_grow2 = strtod(argv[6], NULL);
  expgen2 = strtol(argv[7], NULL, 10);
  sampsize = strtol(argv[8], NULL, 10);
  assert(expgen1 >= 0);
  assert(expgen2 >= 0);
  mut_rate = strtod(argv[9], NULL);
  sel_start = strtol(argv[10], NULL, 10);
  when_bneck = strtol(argv[12], NULL, 10);
  max_site = 40 * startsize * mut_rate * log(startsize);  // will be increased later in program if needed
  nvar = 0;
  count = malloc(max_site * sizeof(int));
  birthgen = malloc(max_site * sizeof(int));
  seed = seed_rng();
  nburn = 10 * startsize;
  //  nburn = startsize / 5;
  //  nburn = 20 * startsize;
  //  nburn = 0;
  fprintf(logf, "max number sites: %d\n", max_site);
  if (sel_start < 0) {sel = sel_store;}
  else {sel = 0;}

  popsize = startsize;
  for (gen = -nburn; gen < 0; gen++) {
    if (gen%1000 == 0) {
      fprintf(logf, "%d, pop: %d, n variants: %d\n", gen, popsize, nvar);
      fprintf(stderr, "%d, pop: %d, n variants: %d\n", gen, popsize, nvar);
    }
    //    next_gen(0.001, sel, gen);    // hack for expanding initial population
    next_gen(0.0, sel, gen);
  }  // end burn-in
  ntot_burn = nvar;
  if (bottlesize > 0 && when_bneck == 0) {
    for (ivar = 0; ivar < nvar; ivar++) {
      freq = (double) count[ivar] / popsize / 2.;
      count[ivar] = brano(freq, 2*bottlesize);
      if (count[ivar] == 0 || count[ivar] == 2*bottlesize) {
	count[ivar] = count[nvar-1];
	birthgen[ivar] = birthgen[nvar-1];
	nvar--;
	ivar--;
      }
    }
    popsize = bottlesize;
    ntot_bneck = nvar;
  }
  r_grow = r_grow1;
  for (gen = 0; gen < 1000000; gen++) {
    if (gen%100 == 0 || (gen > -2 && gen < 10)) {
      fprintf(logf, "%d, pop: %d, n variants: %d\n", gen, popsize, nvar);
      fprintf(stderr, "%d, pop: %d, n variants: %d\n", gen, popsize, nvar);
    }
    if (gen >= sel_start) {sel = sel_store;}
    if (gen == expgen1) {
      r_grow = r_grow2;
      if (bottlesize > 0 && when_bneck == 1) {
	for (ivar = 0; ivar < nvar; ivar++) {
	  freq = (double) count[ivar] / popsize / 2.;
	  count[ivar] = brano(freq, 2*bottlesize);
	  if (count[ivar] == 0 || count[ivar] == 2*bottlesize) {
	    count[ivar] = count[nvar-1];
	    birthgen[ivar] = birthgen[nvar-1];
	    nvar--;
	    ivar--;
	  }
	}
	popsize = bottlesize;
	ntot_bneck = nvar;
      }
    }
    if (gen >= expgen1 + expgen2) {break;}
    next_gen(r_grow, sel, gen);
  }
  if (gen == 1000000) {fprintf(logf, "Exceeded maximum number of generations\n");}
  
  // Draw sample from population and output frequencies
  count_by_class = malloc(2 * popsize * sizeof(int));
  pre_by_class = malloc(2 * popsize * sizeof(int));
  for (ivar = 0; ivar < 2*popsize; ivar++) {
    count_by_class[ivar] = 0;
    pre_by_class[ivar] = 0;
  }
  nvarold = 0;
  for (ivar = 0; ivar < nvar; ivar++) {
    if (sampsize > 0) {
      freq = (double) count[ivar] / popsize / 2.;
      sampled = brano(freq, 2*sampsize);
    }
    else {
      sampled = count[ivar];
    }
    if (sampled == 0 || sampled == 2*popsize) {continue;}
    count_by_class[sampled]++;
    if (birthgen[ivar] < 0) {pre_by_class[sampled]++;}
  }
  strcpy(filename, argv[11]);
  strcat(filename, ".txt");
  outf = fopen(filename, "w");
  if (outf == NULL) {fprintf(logf, "Could not open file %s\n", filename); abort();}
  
  fprintf(logf, "s\t%.5e\nstart size\t%d\n", sel, startsize);
  fprintf(logf, "N gen expansion 1\t%d\nrgrow1\t%.5e\n", expgen1, r_grow1);
  fprintf(logf, "N gen expansion 2\t%d\nrgrow2\t%.5e\n", expgen2, r_grow2);
  fprintf(logf, "final size\t%d\nbottleneck size\t%d\nsample size\t%d\n", popsize, bottlesize, sampsize);
  fprintf(logf, "number of mutations/genome/gen: %.2e\n", mut_rate);
  fprintf(logf, "Final N variants\t%d\n", nvar);
  fprintf(logf, "N variants at end of burn-in\t%d\n", ntot_burn);
  fprintf(logf, "N variants at bottleneck\t%d\n", ntot_bneck);

  // Output absolute frequency spectrum
  for (dafbin = 0; dafbin < ndaf; dafbin++) {
    raw_distr[dafbin] = 0;
    raw_pre[dafbin] = 0;
  }
  pi = 0;
  size = popsize;
  if (sampsize > 0) {size = sampsize;}
  fprintf(stderr, "sampsize: %d\n", size);
  for (iclass = 1; iclass < 2*size; iclass++) {
    if (count_by_class[iclass] == 0) {continue;}
    freq = (double) iclass / size / 2.;
    dafbin = freq / daf_binsize;
    if (dafbin == ndaf) {dafbin = ndaf - 1;}
    raw_distr[dafbin] += count_by_class[iclass];
    raw_pre[dafbin] += pre_by_class[iclass];
    pi += 2 * count_by_class[iclass]* size * 2 * freq * (1-freq) / (2 * size - 1);
  }

  fprintf(outf, "bin\tderived_freq\traw_counts\tpre_bneck");
  fprintf(outf, "\n");
  for (dafbin = 0; dafbin < ndaf; dafbin++) {
    fprintf(outf, "%d\t%.3f\t%d\t%d", dafbin, 
	    (dafbin + 0.5) * daf_binsize, 
	    raw_distr[dafbin], raw_pre[dafbin]);
    fprintf(outf, "\n");
  }

  fprintf(logf, "mut\t%lu\n", totmut);
  fprintf(logf, "pi\t%.4e\n", pi);

  return 0;
}

void next_gen(double r_grow, double sel, int gen) {
  int ivar, nmut, newsize;
  double freq;

  newsize = (int) ((1.+r_grow) * popsize + 0.5);
  // Update old mutations
  for (ivar = 0; ivar < nvar; ivar++) {
    freq = (double) count[ivar] / popsize / 2.;
    if (popsize > 50 && count[ivar] <= 10) {
      count[ivar] = prano( (1-sel)/(1-sel*freq) * count[ivar] * (1+r_grow) );
    }
    else {
      count[ivar] = brano((1-sel)/(1-sel*freq) * freq,  2 * newsize );
    }
    if (count[ivar] == 0 || count[ivar] == 2*newsize) {
      count[ivar] = count[nvar-1];
      birthgen[ivar] = birthgen[nvar-1];
      nvar--;
      ivar--;
    }
  }  // end updating old mutations
  // New mutations
  popsize = newsize;
  nmut = prano(2 * popsize * mut_rate);
  totmut += nmut;

  for (ivar = 0; ivar < nmut; ivar++) {
    count[nvar] = 1;
    birthgen[nvar] = gen;
    nvar++;
    if (nvar == max_site) {
      max_site *= 4;
      if (max_site > 2000000000) {max_site = 2000000000;}
      fprintf(stdout, "Increasing number of sites to %d\n", max_site);
      count = realloc(count, max_site * sizeof(int));
      birthgen = realloc(birthgen, max_site * sizeof(int));
      assert(count != NULL);
      assert(birthgen != NULL);
      if (nvar >= max_site) {
	fprintf(stderr, "Ran out of room in 4-byte int\n");
	abort();
      }
    }
  }  // end new mutations
}
