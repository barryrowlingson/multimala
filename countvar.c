#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define SQR(a) ((a)*(a))
#define ABS(a) ((a) > 0  ?  (a)  : (-a))
#define MIN(a,b) ((a) < (b) ?   (a)  : (b))
#define M 32
#define N 64
#define L1 77.5
#define L2 38.75
#define SIGMA2 1.0
#define BETA 7.5
#define B 0.7
#define MU -9.732664

double rho2(x,y,h)
     double x[2],y[2],h;
{
  return(exp(SIGMA2*exp(-sqrt(SQR(x[0]-y[0])+SQR(x[1]-y[1]))/BETA)
             *exp(-B*h)));
}

void main()
{
  int h,i,j,k,l;
  float tmp;
  double sum,x[2],y[2],lambda[32*64],int1,int2,r,s;
  FILE *infile,*outfile;

  infile=fopen("smoothcount.dat","r");
  for (i=0; i<32; i++)
    for (j=0; j<64; j++) {
      fscanf(infile,"%f",&tmp);
      lambda[i*64+j]=(double)tmp*exp(MU+SIGMA2/2.0)/(L1*L2)*64.0*32.0;
    }
  fclose(infile);
  outfile=fopen("countvar.out","w");
  fprintf(outfile,"Integrated rho^2 over observation window\n");
  fprintf(outfile,"sigma2=%4.2f  beta=%4.2f\n\n",SIGMA2,BETA);
  fprintf(outfile,"   h   Integral\n");
  for (h=0; h<11; h++) {
    sum=0.0;
    for (i=0; i<M; i++)
      for (j=0; j<N; j++)
        for (k=0; k<M; k++)
          for (l=0; l<N; l++)
            {
              x[0]=((double)i+0.5)*L2/(double)M;
              x[1]=((double)j+0.5)*L1/(double)N;
              y[0]=((double)k+0.5)*L2/(double)M;
              y[1]=((double)l+0.5)*L1/(double)N;
              r=32.0/M; s=64.0/N;
              int1=lambda[(int)(i*r)*64+(int)(j*s)];
              int2=lambda[(int)(k*r)*64+(int)(l*s)];
              sum+=rho2(x,y,(double)h)*int1*int2;
            }
    sum=sum*L1*L2*L1*L2/((double)M*(double)N*(double)M*(double)N);
    printf("Sum=%8.4f\n",sum);
    fprintf(outfile,"%4d   %8.6f\n",h,sum);
  }
  fclose(outfile);
}

