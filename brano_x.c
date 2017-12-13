#include <stdio.h>
#include <math.h>
#include "random.h"
#define PI 3.14159265358979
/* binomial deviates, Numerical Recipes */

/* ntrials, prob pp */
int brano(double pp, int n) {  
  double gammln(double);
  int j, bnl;
  static int nold = -1;
  double am, em, g, angle, p, sq, t, y;
  static double pold = -1., pc, plog, pclog, en, oldg;

  if (n == 0) {return 0;}
  p = (pp <= 0.5 ? pp : 1.0 - pp);
  am = n * p;
  if (n < 250) {
    bnl = 0.;
    for (j = 1; j <= n; j++) {
      if (random_double() < p) {
	bnl++;
      }
    }
  }
  else if (am < 1.) {
    g = exp(-am);
    t = 1.;
    for (j = 0; j <= n; j++) {
      t *= random_double();
      if (t < g) {break;}
    }
    bnl = (j <= n) ? j : n;
  }
  else {
    if (n != nold) {
      en = n;
      oldg = gammln(en + 1.);
      nold = n;
    }
    if (p != pold) {
      pc = 1. - p;
      plog = log(p);
      pclog = log(pc);
      pold = p;
    }
    sq = sqrt(2. * am * pc);
    do {
      do {
	angle = PI * random_double();
	y = tan(angle);
	em = sq * y + am;
      } while (em < 0. || em >= (en + 1.));
      em = floor(em);
      t = 1.2 * sq * (1. + y * y) * exp(oldg - gammln(em + 1.) - 
					gammln(en - em + 1.) + em * plog + (en - em) * pclog);
    } while (random_double() > t);
    bnl = (int) (em + 0.5);
  }
  if (p != pp) {
    bnl = n - bnl;
  }
  return bnl;
}

