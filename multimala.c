/*compare this version with Anders', i adjust the seasonal effect*/
  /*and the temporal model is exp(-day/b)*/ 
/* Program for estimation of Y_1,...,Y_k from X_1,...,X_k, where 
   X_1,...,X_k are observations from a space time LGCP with
   OU intensity Y. This version allows an external "population density",
   i.e. the intensity of X_t is exp(Y_t)*lambda, where lambda is a
   user supplied surface. */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define pi 3.141592653589793238462643
typedef float *Field;

void addGauss(Field surf, int dim2, int j, float scale);
void fourn(float data[], unsigned long nn[], int ndim, int isign);
double drand48();

void main(int argc, char **argv){
  int i,j,k,l,nspp,accept,itr,first,N,M,ix,iy,test,dim2,*npts;
  unsigned long nn[2];
  float mu,sigma2,b,*alpha,*beta,*Sigma,*season,*logseason,scale,lx,ly,*sumy,*ssy,*justy,*exceedc2,*exceedc4,*exceedc8,point[3],logfz,logfzprop,qzzprop,qzpropz,divdim2,tmp;
  Field S,y,z,gradz,gradzprop,zprop,yprop,expyprop,u,lambda,loglambda,cellcount,expy,temp;
  FILE *infile1,*infile2,*infile3,*outfile1,*outfile2,*outfile3,*acceptrates,*timeseries;
  char *Sinfn="S.in";
  char *inputfn="multimala.in";
  char *fnaccept="multi.rates.out";
  char *fntime="timeseries.out";
  char gammafn[50],windowfn[50],lambdafn[50],sumyfn[50],ssyfn[50],justfn[50],exceedc2fn[50],exceedc4fn[50],exceedc8fn[50];

 

  /* Print usage if necessary and count the number of spp files. */
  if (argc<5) {
    printf("Usage: multimala sppfile_1 ... sppfile_k first scale itr\n\n");
    exit(1);
  }
  nspp=argc-4;
  alpha=(float *)calloc(nspp,sizeof(float));
  beta=(float *)calloc(nspp,sizeof(float));
  Sigma=(float *)calloc(nspp,sizeof(float));
  season=(float *)calloc(nspp,sizeof(float));
  logseason=(float *)calloc(nspp,sizeof(float));
  npts=(int *)calloc(nspp,sizeof(int));

 
  /* Read input from multimala.in */
  infile1=fopen(inputfn,"r");
  fscanf(infile1,"%f%f%d%d%f%f",&mu,&sigma2,&M,&N,&lx,&ly);
  fscanf(infile1,"%s%s",windowfn,lambdafn);
  if (nspp>1) {
    fscanf(infile1,"%f",&b); /* Check: Do I read in several b's or only one???? */

     for (i=1; i<nspp; i++) {
      fscanf(infile1,"%f",&tmp);
      alpha[i]=1-exp(-tmp/b);
      beta[i]=exp(-tmp/b);
      Sigma[i]=1-exp(-2*tmp/b);
    }
    fclose(infile1);
 }
  M=2*M; N=2*N;
  alpha[0]=1.0; beta[0]=0.0; Sigma[0]=1.0;


  /* Allocation of memory */
  y=(float *)calloc(nspp*M*N*2,sizeof(float));
  z=(float *)calloc(nspp*M*N*2,sizeof(float));
  gradz=(float *)calloc(nspp*M*N*2,sizeof(float));
  gradzprop=(float *)calloc(nspp*M*N*2,sizeof(float));
  zprop=(float *)calloc(nspp*M*N*2,sizeof(float));
  yprop=(float *)calloc(nspp*M*N*2,sizeof(float)); 
  expyprop= (float *)calloc(nspp*M*N*2,sizeof(float));
  expy=(float *)calloc(nspp*M*N*2,sizeof(float));
  cellcount=(float *)calloc(nspp*M*N*2,sizeof(float));
  u=(float *)calloc(M*N*2,sizeof(float));
  lambda=(float *)calloc(M*N*2,sizeof(float));
  loglambda=(float *)calloc(M*N*2,sizeof(float));
  S=(float *)calloc(M*N*2,sizeof(float));
  sumy=(float *)calloc(nspp*M*N,sizeof(float));
  ssy=(float *)calloc(nspp*M*N,sizeof(float));
  justy=(float *)calloc(nspp*M*N,sizeof(float));
  exceedc2=(float *)calloc(nspp*M*N,sizeof(float));
  exceedc4=(float *)calloc(nspp*M*N,sizeof(float));
  exceedc8=(float *)calloc(nspp*M*N,sizeof(float));



  /* Initialisation of variables */
  first=atoi(argv[nspp+1]); scale=atof(argv[nspp+2]); 
  itr=atoi(argv[nspp+3]); 
 
  dim2=M*N; nn[0]=N; nn[1]=M;
  divdim2=dim2; /* this integer to float conversion actually speeds up
		   the program up when divdim2 is used instead of dim2
		   for the divisions below */  
  accept=0; logfz=0.0;
 
    infile1=fopen(Sinfn,"r");
    infile2=fopen(windowfn,"r");
    infile3=fopen(lambdafn,"r");

 
  for (i=0; i<dim2; i++){
 
    fscanf(infile1,"%f",&(S[2*i]));
    fscanf(infile2,"%f",&(u[2*i]));
    fscanf(infile3,"%f",&(lambda[2*i]));

    S[2*i]=sqrt(S[2*i]); 
    u[2*i]*=4.0*lx*ly/(((float)M)*((float)N));
 
   if (lambda[2*i]>0)
      loglambda[2*i]=log(lambda[2*i]);
      else
	loglambda[2*i]=0;

    S[2*i+1]=u[2*i+1]=lambda[2*i+1]=loglambda[2*i+1]=0.0;
       }
    fclose(infile1);
    fclose(infile2);

  for (l=0; l<nspp; l++)
    for (i=0; i<dim2; i++) {
      z[l*M*N*2+2*i]=y[l*M*N*2+2*i]=cellcount[l*M*N*2+2*i]=
	expy[l*M*N*2+2*i]=expyprop[l*M*N*2+2*i]=
	gradzprop[l*M*N*2+2*i]=zprop[l*M*N*2+2*i]=
	gradz[l*M*N*2+2*i]=z[l*M*N*2+2*i+1]=y[l*M*N*2+2*i+1]=
	cellcount[l*M*N*2+2*i+1]=expy[l*M*N*2+2*i+1]=
	expyprop[l*M*N*2+2*i+1]=gradzprop[l*M*N*2+2*i+1]=
	zprop[l*M*N*2+2*i+1]=gradz[l*M*N*2+2*i+1]=0.0;
      sumy[l*M*N+i]=ssy[l*M*N+i]=justy[l*M*N+i]=exceedc2[l*M*N+i]=exceedc4[l*M*N+i]=exceedc8[l*M*N+i]=0.0;
    }







  /* Read point pattern and convert it to a count matrix
     of the same dimensions as the GRF */
  for (l=0; l<nspp; l++) {
    infile1=fopen(argv[l+1],"r");
    test=fscanf(infile1,"%f%f%f",&point[0],&point[1],&point[2]);
    ix=N/2*point[0]/lx; iy=M/2*point[1]/ly;
    cellcount[l*M*N*2+2*(iy*N+ix)]++;
    season[l]=point[2];
    logseason[l]=log(season[l]);
    for (npts[l]=0; test!=EOF; npts[l]++) {
      test=fscanf(infile1,"%f%f%f",&point[0],&point[1],&point[2]);
      ix=N/2*point[0]/lx; iy=M/2*point[1]/ly;
      if (test!=EOF)
	cellcount[l*M*N*2+2*(iy*N+ix)]++;
    }
    fclose(infile1);
    printf("Day %d Seasonal adjust=%f \n",l,season[l]);
  }

  /* Print info */
  if (first)
    printf("\nFirst run. %d iterations.\n",itr);
  else {
    printf("\nUsing results from previous run. %d iterations.\n",itr);
    for (l=0; l<nspp; l++) {
      sprintf(gammafn,"gamma%d.out",l);
      infile1=fopen(gammafn,"r");
      for (i=0; i<dim2; i++) {
	fscanf(infile1,"%f",&(z[l*dim2*2+2*i]));
	y[l*dim2*2+2*i]=z[l*dim2*2+2*i];
      }
      fclose(infile1);
    }
  }
  if (nspp>1)
    printf("mu=%10.6f   sigma^2=%10.6f   scale=%10.4e   b=%10.6f\n",
	   mu,sigma2,scale,b);
  else
    printf("mu=%10.6f   sigma^2=%10.6f   scale=%10.4e\n",
	   mu,sigma2,scale);
  printf("Grid size of extended grid: %dx%d\n",M,N);
  printf("Simulation window: [%6.2f,%6.2f]x[%6.2f,%6.2f]\n",0.0,lx,0.0,ly);
  printf("Observation window read from: %s\n",windowfn);
  printf("%d point processes:\n",nspp);
  for (l=0; l<nspp; l++)
    printf("  %d points read from: %s\n",npts[l],argv[l+1]);
  printf("\n");

  /* Initialise the Markov chain */
  srand48(time(NULL));
  for (l=0; l<nspp; l++) {
    /* y=z and is now transformed back to y */
    fourn(&y[l*dim2*2]-1,nn-1,2,1);
    for (i=0; i<dim2; i++) {
      y[l*dim2*2+2*i]*=S[2*i];
      y[l*dim2*2+2*i+1]*=S[2*i];
    }
    fourn(&y[l*dim2*2]-1,nn-1,2,-1);
    if (l==0)
      for (i=0; i<dim2; i++)
	y[2*i]=sqrt(sigma2)*y[2*i]/divdim2+mu;
    else
      for (i=0; i<dim2; i++)
	y[l*dim2*2+2*i]=sqrt(sigma2)*y[l*dim2*2+2*i]/divdim2+
	  mu*alpha[l]+beta[l]*y[(l-1)*dim2*2+2*i];
  }
  /* Calculate gradient - first step */
  l=nspp-1;
  for (i=0; i<=M/2; i++)
    for (j=0; j<=N/2; j++)
      expy[l*dim2*2+i*N*2+2*j]=cellcount[l*dim2*2+i*N*2+2*j]-
	u[i*N*2+2*j]*season[l]*lambda[i*N*2+2*j]*exp(y[l*dim2*2+i*N*2+2*j]);
  for (l=nspp-2; l>=0; l--)
    for (i=0; i<=M/2; i++)
      for (j=0; j<=N/2; j++)
	expy[l*dim2*2+i*N*2+2*j]=expy[(l+1)*dim2*2+i*N*2+2*j]*beta[l+1]+
	  cellcount[l*dim2*2+i*N*2+2*j]-u[i*N*2+2*j]*season[l]*lambda[i*N*2+2*j]*
	  exp(y[l*dim2*2+i*N*2+2*j]);
  for (l=0; l<nspp; l++) {
    /* Calculate gradient - second step */
    fourn(&expy[l*dim2*2]-1,nn-1,2,1);
    for (i=0; i<dim2; i++) {
      expy[l*dim2*2+2*i]*=S[2*i];
      expy[l*dim2*2+2*i+1]*=S[2*i];
    }
    fourn(&expy[l*dim2*2]-1,nn-1,2,-1);
    for (i=0; i<dim2; i++)
      gradz[l*dim2*2+2*i]=sqrt(sigma2)*expy[l*dim2*2+2*i]/divdim2-
	z[l*dim2*2+2*i]/Sigma[l];
    /* Calculate log density y */
    for (i=0; i<dim2; i++)
      logfz-=0.5*pow(z[l*dim2*2+i*2],2)/Sigma[l];
    for (i=0; i<=M/2; i++)
      for (j=0; j<=N/2; j++)
	logfz+=cellcount[l*dim2*2+i*N*2+2*j]
	  *(y[l*dim2*2+i*N*2+2*j]+loglambda[i*N*2+2*j]+logseason[l])-
		u[i*N*2+2*j]*lambda[i*N*2+2*j]*season[l]*exp(y[l*dim2*2+i*N*2+2*j]);
  }

  /* The Markov chain iterations begin */
  acceptrates=fopen(fnaccept,"w");
  timeseries=fopen(fntime,"w");
  for (k=0; k<itr; k++){
    for (l=0; l<nspp; l++) {
      /* Simulate proposal zprop */
      for (i=0; i<dim2/2; i++){
	zprop[l*dim2*2+i*2]=z[l*dim2*2+i*2]+0.5*gradz[l*dim2*2+2*i]*scale;
	zprop[l*dim2*2+dim2+i*2]=z[l*dim2*2+dim2+i*2]+
	  0.5*gradz[l*dim2*2+dim2+2*i]*scale;
	addGauss(&zprop[l*dim2*2],dim2,i,scale);
      }
      /* Generate yprop from zprop */
      for (i=0; i<dim2; i++) {
	yprop[l*dim2*2+2*i]=zprop[l*dim2*2+2*i];
	yprop[l*dim2*2+2*i+1]=0.0;
      }
      fourn(&yprop[l*dim2*2]-1,nn-1,2,1);
      for (i=0; i<dim2; i++) {
	yprop[l*dim2*2+2*i]*=S[2*i];
	yprop[l*dim2*2+2*i+1]*=S[2*i]; 
      }
      fourn(&yprop[l*dim2*2]-1,nn-1,2,-1);
      if (l==0)
	for (i=0; i<dim2; i++){
	  yprop[2*i]=sqrt(sigma2)*yprop[2*i]/divdim2+mu;
	  expyprop[2*i]=0.0; expyprop[2*i+1]=0.0;
	}
      else
	for (i=0; i<dim2; i++){
	  yprop[l*dim2*2+2*i]=sqrt(sigma2)*yprop[l*dim2*2+2*i]/divdim2+
	    mu*alpha[l]+beta[l]*yprop[(l-1)*dim2*2+2*i];
	  expyprop[l*dim2*2+2*i]=0.0; expyprop[l*dim2*2+2*i+1]=0.0;
	}
    }
    /* Calculate gradient - first step */
    l=nspp-1;
    for (i=0; i<=M/2; i++)
      for (j=0; j<=N/2; j++)
	expyprop[l*dim2*2+i*N*2+2*j]=cellcount[l*dim2*2+i*N*2+2*j]-
	  u[i*N*2+2*j]*lambda[i*N*2+2*j]*season[l]*exp(yprop[l*dim2*2+i*N*2+2*j]);
    for (l=nspp-2; l>=0; l--)
      for (i=0; i<=M/2; i++)
	for (j=0; j<=N/2; j++)
	  expyprop[l*dim2*2+i*N*2+2*j]=
	    expyprop[(l+1)*dim2*2+i*N*2+2*j]*beta[l+1]+
	    cellcount[l*dim2*2+i*N*2+2*j]-u[i*N*2+2*j]*lambda[i*N*2+2*j]*
	        season[l]*exp(yprop[l*dim2*2+i*N*2+2*j]);
    logfzprop=0.0; qzzprop=0.0; qzpropz=0.0;
    for (l=0; l<nspp; l++) {
      /* Calculate gradient - second step */
      fourn(&expyprop[l*dim2*2]-1,nn-1,2,1);
      for (i=0; i<dim2; i++){
	expyprop[l*dim2*2+2*i]*=S[2*i]; 
	expyprop[l*dim2*2+2*i+1]*=S[2*i]; 
      }
      fourn(&expyprop[l*dim2*2]-1,nn-1,2,-1);
      for (i=0; i<dim2; i++)
	gradzprop[l*dim2*2+2*i]=sqrt(sigma2)*expyprop[l*dim2*2+2*i]/divdim2-
	  zprop[l*dim2*2+2*i]/Sigma[l];
      /* Calculate log density for proposal */
      for (i=0; i<dim2; i++)
	logfzprop-=0.5*pow(zprop[l*dim2*2+i*2],2)/Sigma[l];
      for (i=0; i<=M/2; i++)
	for (j=0; j<=N/2; j++)
	  logfzprop+=cellcount[l*dim2*2+i*N*2+2*j]
	    *(yprop[l*dim2*2+i*N*2+2*j]+loglambda[i*N*2+2*j]+logseason[l])
	    -u[i*N*2+2*j]*lambda[i*N*2+2*j]*season[l]*exp(yprop[l*dim2*2+i*N*2+2*j]);
      /* Calculate acceptance probabilities */
      for (i=0; i<dim2; i++){
	qzzprop+=pow(zprop[l*dim2*2+2*i]-
		     (z[l*dim2*2+2*i]+0.5*gradz[l*dim2*2+2*i]*scale),2);
	qzpropz+=pow(z[l*dim2*2+2*i]-
		     (zprop[l*dim2*2+2*i]+
		      0.5*gradzprop[l*dim2*2+2*i]*scale),2);
      }
    }
    qzzprop*=(-0.5/scale); qzpropz*=(-0.5/scale);
    /* Check whether to accept proposal or not */
    if (log(drand48())<logfzprop+qzpropz-logfz-qzzprop){
      logfz=logfzprop;
      temp=gradz;
      gradz=gradzprop;
      gradzprop=temp;
      temp=y;
      y=yprop;
      yprop=temp;
      temp=z;
      z=zprop;
      zprop=temp;
      accept++;
    }
    /* Update sample sums and ss's */
    for (l=0; l<nspp; l++)
      for (i=0; i<dim2; i++){
	sumy[l*dim2+i]+=y[l*dim2*2+2*i];
	ssy[l*dim2+i]+=pow(y[l*dim2*2+2*i],2);
	justy[l*dim2+i]=y[l*dim2*2+2*i];
          if (y[l*dim2*2+2*i]>log(2))
               {exceedc2[l*dim2+i]=exceedc2[l*dim2+i]+1;}
          if (y[l*dim2*2+2*i]>log(4))
               {exceedc4[l*dim2+i]=exceedc4[l*dim2+i]+1;}
          if (y[l*dim2*2+2*i]>log(8))
               {exceedc8[l*dim2+i]=exceedc8[l*dim2+i]+1;}         
      }


    /* Save the value of selected pixels every 10 iterations */
    if ((k+1)%10==0){
      fprintf(timeseries," % 15.10e % 15.10e % 15.10e % 15.10e % 15.10e\n",
	      y[26*N*2+2*10],y[9*N*2+2*18],y[13*N*2+2*47],
	      y[30*N*2+2*51],y[26*N*2+2*34]);
      fflush(timeseries);
    }
    /* Print the acceptance rate every 100 iterations */
    if ((k+1)%100==0){
      fprintf(acceptrates,"Itr: %d rate: %f\n",k+1,accept/100.0);
      fflush(acceptrates);
      printf("Itr: %d rate: %4.2f\n",k+1,accept/100.0);
      accept=0;
    }
    if ((k+1)%1000==0){
      for (l=(nspp-1); l<nspp; l++){
	sprintf(justfn,"justy%d_%d_%d.out",nspp,l,k+1);
	outfile1=fopen(justfn,"w");
	  for (i=0; i<dim2; i++){
	    fprintf(outfile1," % 15.10e ",justy[l*dim2+i]);
	  }
	fclose(outfile1);
      }
    }
  }
  fclose(acceptrates); fclose(timeseries);

  /* Write results to files */
  for (l=0; l<nspp; l++) {
    sprintf(gammafn,"gamma%d.out",l);
    sprintf(sumyfn,"sumy%d_%d.out",nspp,l);
    sprintf(ssyfn,"ssy%d_%d.out",nspp,l);
    outfile1=fopen(gammafn,"w");
    outfile2=fopen(sumyfn,"w");
    outfile3=fopen(ssyfn,"w");
    for (i=0; i<dim2; i++){
      fprintf(outfile1," % 15.10e ",z[l*dim2*2+2*i]);
      fprintf(outfile2," % 15.10e ",sumy[l*dim2+i]);
      fprintf(outfile3," % 15.10e ",ssy[l*dim2+i]);
    }
    fclose(outfile1); fclose(outfile2); fclose(outfile3);
  }
  for (l=(nspp-1); l<nspp; l++){
     sprintf(exceedc2fn,"exceed%d_%d_2.out",nspp,l);
     sprintf(exceedc4fn,"exceed%d_%d_4.out",nspp,l);
     sprintf(exceedc8fn,"exceed%d_%d_8.out",nspp,l);
     outfile1=fopen(exceedc2fn,"w");
     outfile2=fopen(exceedc4fn,"w");
     outfile3=fopen(exceedc8fn,"w");
     for (i=0; i<dim2; i++){
        fprintf(outfile1," % 15.10e ",exceedc2[l*dim2+i]);
        fprintf(outfile2," % 15.10e ",exceedc4[l*dim2+i]);
        fprintf(outfile3," % 15.10e ",exceedc8[l*dim2+i]);
     }
     fclose(outfile1);fclose(outfile2); fclose(outfile3);
  }

}



void addGauss(Field surf, int dim2, int i, float scale)
     /* Add an N(0,scale) distributed error to elements 2*i 
	and 2*i+dim2 of surf */
{
  double theta,E,R;
  
  theta=2*pi*drand48();
  E=-log(drand48());
  R=sqrt(2*E);
  surf[2*i]+=sqrt(scale)*R*sin(theta);
  surf[2*i+dim2]+=sqrt(scale)*R*cos(theta);
}

void fourn(float data[], unsigned long nn[], int ndim, int isign)
     /* Fast Fourier transfrom from Numerical Recipes for C */
{
  int idim;
  unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
  unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
  float tempi,tempr;
  double theta,wi,wpi,wpr,wr,wtemp;
  
  for (ntot=1,idim=1;idim<=ndim;idim++)
    ntot *= nn[idim];
  nprev=1;
  for (idim=ndim;idim>=1;idim--) {
    n=nn[idim];
    nrem=ntot/(n*nprev);
    ip1=nprev << 1;
    ip2=ip1*n;
    ip3=ip2*nrem;
    i2rev=1;
    for (i2=1;i2<=ip2;i2+=ip1) {
      if (i2 < i2rev) {
	for (i1=i2;i1<=i2+ip1-2;i1+=2) {
	  for (i3=i1;i3<=ip3;i3+=ip2) {
	    i3rev=i2rev+i3-i2;
	    SWAP(data[i3],data[i3rev]);
	    SWAP(data[i3+1],data[i3rev+1]);
	  }
	}
      }
      ibit=ip2 >> 1;
      while (ibit >= ip1 && i2rev > ibit) {
	i2rev -= ibit;
	ibit >>= 1;
      }
      i2rev += ibit;
    }
    ifp1=ip1;
    while (ifp1 < ip2) {
      ifp2=ifp1 << 1;
      theta=isign*6.28318530717959/(ifp2/ip1);
      wtemp=sin(0.5*theta);
      wpr = -2.0*wtemp*wtemp;
      wpi=sin(theta);
      wr=1.0;
      wi=0.0;
      for (i3=1;i3<=ifp1;i3+=ip1) {
	for (i1=i3;i1<=i3+ip1-2;i1+=2) {
	  for (i2=i1;i2<=ip3;i2+=ifp2) {
	    k1=i2;
	    k2=k1+ifp1;
	    tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
	    tempi=(float)wr*data[k2+1]+(float)wi*data[k2];
	    data[k2]=data[k1]-tempr;
	    data[k2+1]=data[k1+1]-tempi;
	    data[k1] += tempr;
	    data[k1+1] += tempi;
	  }
	}
	wr=(wtemp=wr)*wpr-wi*wpi+wr;
	wi=wi*wpr+wtemp*wpi+wi;
      }
      ifp1=ifp2;
    }
    nprev *= n;
  }
}

