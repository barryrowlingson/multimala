/* Simulates an inhomogeneous Poisson process in a rectangular window. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "ranlib.h" /*In Anders' code, <ranlib.h>*/

void sim_pois();

void sim_pois(double *intensity,
	      long *npts,
	      double *l1,
	      double *l2,
	      long *N,
	      long *M,
	      double *spp)
{
  int i,j,k,l;
  double maxint,u,pitch[2];
  
  /*
  printf("Simulating inhomogeneous Poisson process with %d points.\n",*npts);
  */
  pitch[0]=*l1/(*M);
  pitch[1]=*l2/(*N);
  maxint=intensity[0];
  /* maxx=pitch[0]/2.0; maxy=pitch[1]/2.0;*/
  for (i=0; i<*N; i++)
    for (j=0; j<*M; j++)
      if (intensity[i*(*M)+j]>maxint) {
	maxint=intensity[i*(*M)+j];
      }

  for (i=0; i<*npts; i++) {
    do {
      spp[2*i]=genunf(0,*l1);
      spp[2*i+1]=genunf(0,*l2);
      u=genunf(0,1);
      k=(int)(floor(spp[2*i]/pitch[0]));
      l=(int)(floor(spp[2*i+1]/pitch[1]));
    } while ((u>=intensity[l*(*M)+k]/maxint) || (k>=128) || (l>=128));

  }
}
