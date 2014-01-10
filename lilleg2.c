
#include<math.h>
#include <stdlib.h>
#include <stdio.h>

#define PI 3.1415927
#define SQR(a) ((a)*(a))
#define ABS(a) ((a) > 0  ?  (a)  : (-a))
#define MIN(a,b) ((a) < (b) ?   (a)  : (b))
#define STEPS 200


/* Prototype*/

double fast_epker(double t,double h);
double dist(double x1,double x2,double y1,double y2);
double k1(double d1,double d2,double dist);
double k2(double d1,double d2,double dist);
double ripley(double x, double y, double d, int *isinw,
	      int m, int n, double l1,double l2);
void ns_lilleg(double *data1, double *data2, double *u,
	       int *m, int *n, double *g, double *eph,
	       double *e1, double *e2, double *tlambda,
	       int *tisinw, int *nx, int *ny, double *area);
void min_dist (double *spp,int *npts,double *mindist);
void covequn (double *lambdadx, long int *m,     long int *n,
	     double *w,  double *v, double *vday,  int *maxvday,
	     double *sigma2,double *beta,double *phi, double *res);
void distxy (double *grid, int *z, double *res);




double fast_epker(double t,double h)
{
  return((0.75/h)*(1-SQR(t/h)));
}


double dist(double x1,double x2,double y1,double y2)
{
  double dd;
  
  dd=sqrt(SQR(x1-x2)+SQR(y1-y2));
  return(dd);
} 


double k1(double d1,double d2,double dist)
{
  double arg1,arg2,kk;
  
  if (dist<=d1)
    arg1=0;
  else
    arg1=acos(d1/dist);
  if (dist<=d2)
    arg2=0;
  else
    arg2=acos(d2/dist);
  kk=1/(1-((arg1+arg2)/PI));
  return(kk);
}


double k2(double d1,double d2,double dist)
{
  double arg1,arg2,kk;

  arg1=acos(d1/dist);
  arg2=acos(d2/dist);
  kk=1/(0.75-(0.5*(arg1+arg2)/PI));
  return(kk);
}


double ripley(double x, double y, double d, int *isinw,
                 int m, int n, double l1,double l2)
{
  int i,k,l;
  double u,v,theta,weight,res;
 

  weight=0.0;
  for (i=0; i<STEPS; i++) {
    theta=(double)((i*2.0*PI)/STEPS);

    u=x+d*cos(theta);
    v=y+d*sin(theta);
    l=(int)(floor(u*n/l1));
    k=(int)(floor(v*m/l2));

    if ((l>=0) && (l<n) && (k>=0) && (k<m)) {
      if (isinw[k*n+l]>0)
           weight++;
    }
  }
  res=(double)(weight/STEPS);
  /*printf("weight: %f \n",res);*/
  return(res);
}




/*data1:x coordinate;  data2:y coordinate*/
/*u:distance of interested;*/ 
/*m: how many distance we are interested*/
/*n:how many data point in this day*/
/*eph:EP bandwith;  e1e2:l1l2 */
/*tlambda and tisinw :transpose the original lambda and isinw to match nrow=length(x), ncol=length(y)*/

/*ix and iy are defined as integer.  rounded down to interger here*/
/*so suppose ix is in grid (81,30), will generate (ix,iy)(80,29)*/
/*this is good since when convert matrix lamda from R to C*/
/*the index starts from 0 to n-1, this is the right index to be used!!!*/
void ns_lilleg(double *data1, double *data2, double *u,
	       int *m, int *n, double *g, double *eph,
	       double *e1, double *e2, double *tlambda,
	       int *tisinw, int *nx, int *ny, double *area)
{

  double t,t1,tmp,abc;
  double rho3,help;
  int i,ix,iy,j,l;
  
for (l=0;l<(*m);l++) {
    rho3=0.0;
    g[l]=0.0;

    for(i=0;i<(*n);i++) {
      for(j=0;j<(*n);j++) { 

        if(j!=i){

  
          t1=dist(data1[i],data1[j],data2[i],data2[j]);
          t=(u[l]-t1);
  
       
	  if (fabs(t)<(*eph)) {

	    ix=(int)(floor(data1[i]*(*nx)/(*e1)));
	    iy=(int)(floor(data2[i]*(*ny)/(*e2)));
	    tmp=tlambda[iy*(*nx)+ix];
           if (tmp<0.00000000000000001) 
                 {printf("intensity shouldn't be zero \n");
		 printf("DATA i:(%f,%f), put on C grid(%d,%d) tmp=%f \n",data1[i],data2[i],ix,iy,tmp);}

	    ix=(int)(floor(data1[j]*(*nx)/(*e1)));
	    iy=(int)(floor(data2[j]*(*ny)/(*e2)));
	    tmp=tmp*tlambda[iy*(*nx)+ix];

            if (tmp<0.00000000000000001) 
                 {printf("intensity is zero \n");
		  printf("DATA j:(%f,%f), put on C grid(%d,%d), tmp=%f \n",data1[j],data2[j],ix,iy,tlambda[iy*(*nx)+ix]);
                  }

       	    help=fast_epker(t,(*eph))/tmp;
	    abc=ripley(data1[i],data2[i],t1,tisinw,(*ny),(*nx),(*e1),(*e2));
	   
	    rho3=rho3+help/abc;

	  } 
	}
      }
    }

    g[l]=rho3/(2*PI*u[l]*(*area));

  }
}



void min_dist(double *spp, int *npts, double *mindist)
{
  int i,j;
  double tmp,a,b;

  a=spp[0]-spp[2];
  b=spp[1]-spp[3];
  *mindist=sqrt(a*a+b*b);
  for (i=0; i<*npts; i++)
    for (j=i+1; j<*npts; j++) {
        a=spp[2*i]-spp[2*j];
	b=spp[2*i+1]-spp[2*j+1];
	tmp=sqrt(a*a+b*b);
	if (tmp<*mindist)
	  *mindist=tmp;
    }
}



void covequn(double *lambdadx, long int *m,     long int *n,
	     double *w,  double *v, double *vday,  int *maxvday,
	     double *sigma2,double *beta,double *phi, double *res)
{  
  double int1,int2,spatsig2,distant,u;
  register long int i,j,k,l,p;

  for (p=0; p<(*maxvday); p++){
    res[p]=0.0;
    spatsig2=(*sigma2)*exp(-(vday[p])/(*beta));
    printf("Day p: %f \n",vday[p]);
    for (i=0; i<*m; i++){
        for (j=0; j<*n; j++){
	 int1=lambdadx[j*(*m)+i];
        for (k=0; k<*m; k++){
	  u=v[i*(*m)+k];
          for (l=0; l<*n; l++)
            {
              int2=lambdadx[l*(*m)+k];
	      distant=sqrt( w[j*(*n)+l] + u );
	      res[p]+=int1*int2*(exp(spatsig2*exp(-distant/(*phi)) ));
	    }
	  }
	
       }
    }
    /* Rprintf("C(V;beta): %f \n",res[p]);*/
  }
}


void distxy (double *grid, int *z, double *res)

{

  int j,l,ind;

   for (j=0; j<*z; j++){
       for (l=0; l<*z; l++){
	 ind=j*(*z)+l;
	 res[ind]=0.0;
	 res[ind]=(grid[j]-grid[l])*(grid[j]-grid[l]);
       }
   }

}






