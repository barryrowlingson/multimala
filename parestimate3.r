library(splancs)
library(shapelib)
library(ts)
library(adapt,lib.loc="~/Rlibs")
library(sn,lib.loc="~/Rlibs")
dyn.load(paste(UTILITYDIR,"lilleg2.so",sep=""))
dyn.load(paste(UTILITYDIR,"r.poisson.sim.so",sep=""))
dyn.load(paste(UTILITYDIR,"libadaptpoly.so",sep=""))


inarea <- function(spp.x,spp.y,e1,e2,nx,ny,tisinw,xgrid,ygrid){
  k <- floor((spp.x-xgrid[1])/e1*128)+1
  l <- floor((spp.y-ygrid[1])/e1*128)+1
  if (k>128 | k <1 | l>128 | l<1)
    {return(0)}
  else{
    return(tisinw[k,l])}
}  

#data must be in 3 column, 1st-x cor in meter; 2nd-y cor; 3rd alg data
sppdata <- function(data,xrange,yrange){
  data$XAxis <-data[,1]-xrange[1]
  data$YAxis <- data[,2]-yrange[1]
  data$ALGDATE_strptime(as.character(data[,3]),format="%d/%m/%Y")
  data$SETNO_floor(julian(data$ALG,origin=as.POSIXct("2000-07-31")))
  data$WDAY2_as.POSIXlt(data$ALG)$wday #sun to sat get 0 to 6
  Fweekdy <- c(2.62774,-0.27164,-0.47594,-0.36660
               ,-0.43936 ,-0.44552,-0.12456,-0.13525
               ,0.04142,-0.15387 ,0.06523)

  W <- 2*pi/365
  data$Time <- Fweekdy[data$WDAY2+1]+
            Fweekdy[8]*cos(W*data$SETNO)+Fweekdy[9]*sin(W*data$SETNO)+
            Fweekdy[10]*cos(2*W*data$SETNO)+Fweekdy[11]*sin(2*W*data$SETNO)
  if(!(sum(data$WDAY2)==0)) {data$Time <- data$Time+Fweekdy[1]}
  data$byglm <- exp( data$Time)
  return(cbind(data$XAxis,data$YAxis,data$byglm))
 # write.table(spp,file="NHSOCT20.txt", sep=" ",row.names=F, col.names=F)
}


#### spatial (two-dimensional) kernel density estimation ##########


kernel.2d <- function(data,x,y,h,kernel="normal",weight,intensity=F){
  z <- matrix(0,nrow=length(x),ncol=length(y))
  if(kernel=="normal"){
    for(i in seq(along=data[,1])){
      z <- z+weight[i]*as.matrix(dnorm((x-data[i,1])/h)/h)%*%
        t(as.matrix(dnorm((y-data[i,2])/h)/h))
    NULL
    }
   } 
  if(kernel=="uniform"){
    hmat <- matrix(h,nrow=length(x),ncol=length(y))
    for(i in seq(along=data[,1])){

      dmat <- sqrt(matrix((x-data[i,1])^2,byrow=F,nrow=length(x),ncol=length(y))+
                matrix((y-data[i,2])^2,byrow=T,nrow=length(x),ncol=length(y)))
      z <- z+weight[i]*as.numeric(dmat<=h)/(pi*h^2)
    NULL
    }
  }
  if(kernel=="quartic"){
   for(i in seq(along=data[,1])){

       distmat <- matrix((x-data[i,1])^2,byrow=F,nrow=length(x),ncol=length(y))+
               matrix((y-data[i,2])^2,byrow=T,nrow=length(x),ncol=length(y))

       kdensmat <- (3/(8*pi)*(1-distmat/(8*(h^2)))^2)/(h^2)
       kdensmat[distmat>(8*h*h)] <- 0
       z <- z+weight[i]*kdensmat                                     
   NULL
     }
   }

 if(kernel=="quadratic"){
   for(i in seq(along=data[,1])){

       distmat <- matrix((x-data[i,1])^2,byrow=F,nrow=length(x),ncol=length(y))+
                matrix((y-data[i,2])^2,byrow=T,nrow=length(x),ncol=length(y))

       kdensmat <- ((2/pi)*(1-distmat/(h^2)))/(h^2)
       kdensmat[distmat>(h*h)] <- 0
       z <- z+weight[i]*kdensmat                                     
   NULL
     }
 }
                                         
if(intensity){ans <- z}
  else{ans <- z/nrow(data)}
ans
}

#this is the function to get each grid point from the planner
#result will be a 2 columns matrix goes by
#(x1,y1)(x2,y1)...(xn,y1)... (x1,yn)(x2,yn)...(xn,yn)
grid.cor<-function(xgrid,ygrid){
  ygrid2 <- NULL

  for (i in 1:length(ygrid)){
    ygrid2<-c(ygrid2,rep(ygrid[i],length(xgrid)))
  }
  ans<- cbind(rep(xgrid,length(ygrid)),ygrid2)
ans
}

#adjust constant for the boundry effect
#x,y grids to calculate kernel

adaptpoly<-function(poly, xpts, h, kernel, c1, c2, eps, mcalls){
     if(is.null(nrow(xpts))){
         ans<-.C("adaptpoly", as.double(poly), as.integer(nrow(poly)),
                 as.double(xpts), as.integer(1), as.double(h), 
                 as.integer(kernel),as.double(c1), as.double(c2),
                 as.double(c(range(poly[,1]),range(poly[,2]))), 
                 as.double(eps),err=as.double(0), as.integer(mcalls),
                 ncalls=as.integer(0),
                 ier=integer(6), cxh=as.double(0))
     } else {
         ans<-.C("adaptpoly", as.double(poly), as.integer(nrow(poly)),
                 as.double(xpts), as.integer(nrow(xpts)), as.double(h), 
                 as.integer(kernel),as.double(c1), as.double(c2),
                 as.double(c(range(poly[,1]),range(poly[,2]))), 
                 as.double(eps),err=double(nrow(xpts)),
                 as.integer(mcalls), ncalls=integer(nrow(xpts)),
                 ier=integer(6), cxh=double(nrow(xpts)))
     }
     list(err=ans$err, ncalls=ans$ncalls, ier=ans$ier, cxh=ans$cxh)

   }


#The muest under the STATIONARY case (with no time trend adjusted),
#is a constant and length=period
#in with-time-adjustment senario, this is a vector with length=period.

gfunc <- function(data,period,minsetno,xpts3,xrange,yrange,EPC,l1,l2,
                    muest,isinw,poplambda){

  ypts <- matrix(0,nrow=period+1,ncol=length(xpts3))
  days <- rep(0,length=period+1)

  aven <- dim(data)[1]/period
  area <- sum(isinw)*l1*l2/(prod(dim(isinw)))

  EPrho <- aven/area
  EPh <- (EPC)/sqrt(EPrho)

  
 for (i in 1:period)   # Here the estimated mean intensity is used
  {

    if (i%%100==0) {cat("done",i,"\n")}

    if (sum( data$SETNO==(minsetno-1+i) )==0)
      { cat("No Case in SETNO",(minsetno-1+i),"\n")
        ypts[i+1,]<- ypts[i,]
        days[i+1] <- days[i]
      }

    else
      {
      
       X<-cbind((data$XAxis[data$SETNO==(minsetno-1+i)]-xrange[1]),
                (data$YAxis[data$SETNO==(minsetno-1+i)]-yrange[1]))

       inten <- muest[i]*poplambda

       ypts[i+1,] <- ypts[i,]+ns.g.est(X,xpts3,EPh,l1,l2,isinw,area,inten)
      days[i+1] <- days[i]+1
       
      }
 }

  ypts2 <- ypts[(period+1),]/days[(period+1)]  
  return(ypts2)

}




#copy from Anders' code but take off the rho since that doesn't use anywhere.
#xpts3 is a range of distance interested
#EPh is Ep bandwith-based on intensity
#l1 is length of x/ l2 is length of y
#isinw is a 0/1 matrix which nrow=length(y), ncol=length(x), indicate study area.
#lambdaest is the intensity at this day, which the coordinate matched isinw.

ns.g.est <- function(spp.data,u,EPh,l1,l2,isinw,area,lambdaest) {
  x <- as.double(spp.data[,1])
  y <- as.double(spp.data[,2])
  n <- as.integer(length(x))
  m <- as.integer(length(u))
  storage.mode(u) <- "double"
  storage.mode(EPh) <- "double"
  storage.mode(l1) <- "double"
  storage.mode(l2) <- "double"
  tisinw <- as.integer(t(isinw))
  tlambda <- as.double(t(lambdaest))
  nx <- as.integer(dim(lambdaest)[2])
  ny <- as.integer(dim(lambdaest)[1])
  storage.mode(area) <- "double"
  res <- .C("ns_lilleg",x,y,u,m,n,g=double(length(u)),EPh,l1,l2,tlambda,tisinw,nx,ny,area)$g

 return(res)

}




#udist is the separation of the distance


equn.14 <- function(parameter,ghat,udist,epsilon,cons){
  temp <- c(epsilon,udist[1:(length(udist)-1)])
  du <- udist-temp
  browser()
  sigma2 <- parameter[1]
  phi <- parameter[2]
  logg <- log(ghat)
  covfn <- exp(-udist/phi)
  sigphi <- sum(  ((logg)^cons - (sigma2*covfn)^cons)^2 *du)
  return(sigphi)
}

#optimization to get phi, sigma2
#cons is the constant we try to stablize the min contrast estimator
#lambda:space intensity

par.est <- function(start,equn.14,ghat,udist,epsilon,cons,pixelarea){
  equn.14.min <- optim(start,equn.14,method="Nelder-Mead",
		ghat=ghat,udist=udist,epsilon=epsilon,
		cons=cons,control=list(maxit=200, trace=2))

  return(list(sigmasqr=equn.14.min$par[1],phispace=equn.14.min$par[2]))
}

#Chat: temporal covariance: w/o time adjustment
#dailycount: 2 columns; 1st daily count/2nd DAY-OF-THE WEEK
#expcount:  average count for day of the week;
#two columns, 1st:day of the week (sun-sat:0-6) 2nd:average 
Chat <- function(dailycount,vday,expcount){
  Chat <- rep(0,length(vday))

  for (j in 1:length(vday)){
    a <- 0
    for (i in 1:(length(dailycount[,1])-vday[j])){
      if (dailycount[i,1] !=0 & dailycount[i+vday[j],1] !=0){
         
        tempC <- dailycount[i,1]*dailycount[i+vday[j],1]-expcount[(dailycount[i,2]+1),2]*expcount[(dailycount[(i+vday[j]),2]+1),2]

        Chat[j] <- Chat[j]+tempC
         a <- a+1
       }
    }    
    Chat[j] <- Chat[j]/a
  }  

  return(list(Chat=Chat))
               }


Chat.t <- function(dailycount,vday,expcount){
  Chat <- matrix(0,nrow=length(dailycount),ncol=length(vday))
  timewt <- matrix(0,nrow=length(dailycount),ncol=length(vday))
  a <- rep(0,length(vday))
  
  for (j in 1:length(vday)){
    for (i in 1:(length(dailycount)-vday[j])){
        
      if (dailycount[i] !=0 & dailycount[i+vday[j]] !=0){

        timewt[i,j] <- expcount[i]*expcount[i+vday[j]]
        
        Chat[i,j] <- dailycount[i]*dailycount[i+vday[j]]-timewt[i,j]
         a[j] <- a[j]+1
       }
     }    
   }
  return(list(Chat=Chat,timewt=timewt,a=a))
}


#optimization of beta
#WHILE INPUT CASES here, input with 2 col. 1st cases, 2nd DOW
#WHILE INPUT WKCASE, 2 COL, FIRST DOW, 2ND MEAN CASES FOR DOW
beta.est <- function(start,vday,sigma2,phi,epsilon,lambda,pixelarea,xgrid,ygrid,Chat,AVECASE){
  #TRY TO GENERATE A DISTANCE VECTOR (X1-X2)^2 used in cov.v.beta TO FASTEN THE OPTIM
  n <- length(xgrid)
  m <- length(ygrid)

  w <- .C("distxy",as.double(xgrid),as.integer(n),res=double(n*n))$res
  v <- .C("distxy",as.double(ygrid),as.integer(m),res=double(m*m))$res

  lambdadx<-lambda*pixelarea
  dv <- rep(1,length(vday))


  newequn.15 <- function(parameter,Chat,lambdadx,w,v,vday,sigma2,phi,dv,m,n,AVECASE){
    beta <- parameter[1]

     tmp <- cov.v.beta(lambdadx,w,v,vday,sigma2,beta,phi,m,n)
     tmp2 <- (AVECASE)^2*(tmp-1)
     tmp2[1] <- tmp2[1]+AVECASE
     tmp.cor <- tmp2/tmp2[1]
    
     print(cbind(beta,tmp.cor,tmp)) 
     return(sum((Chat-tmp.cor)^2*dv))

  }
  
 equn.15.min <- optim(start,newequn.15,method="Nelder-Mead",
                   control=list(maxit=200, trace=2),
                   Chat=Chat,lambdadx=lambdadx,
                   w=w,v=v,vday=vday,sigma2=sigma2,
                   phi=phi,dv=dv,m=m,n=n,
                   AVECASE=AVECASE)

  return(list(beta=equn.15.min$par[1]))
}



##the beta estimate fun.  consider the week of the day effect!
#beta.est.t <- function(start,vday,sigma2,phi,lambda,pixelarea,xgrid,ygrid,expbyglm,Chat,timewt){
#TRY TO GENERATE A DISTANCE VECTOR (X1-X2)^2 used in cov.v.beta TO FASTEN THE OPTIM

beta.est.t <- function(vday,sigma2,phi,lambda,pixelarea,xgrid,ygrid,expbyglm,Chat,timewt,RATIO){
  n <- 128
  m <- 128

  w <- .C("distxy",as.double(xgrid),as.integer(128),res=double(128*128))$res
  v <- .C("distxy",as.double(ygrid),as.integer(128),res=double(128*128))$res
  lambdadx <- as.matrix(lambda*pixelarea)
  
  dv <- rep(1,length(vday))

  print("start beta.est.t")
  newequn.15 <- function(parameter,Chat,lambdadx,w,v,vday,sigma2,phi,dv,m,n,timewt,expbyglm,RATIO){
    beta <- parameter[1]
    cat("beta=",beta,"\n")
    tmp1 <- matrix(0,nrow=dim(Chat)[1],ncol=dim(Chat)[2])
    tmp <- cov.v.beta(lambdadx,w,v,vday,sigma2,beta,phi,m,n)
    print(cbind(beta,vday,tmp))
    
    for (j in 1:length(vday)){
      tmp1[,j] <- timewt[,j]*(tmp[j]-1)*RATIO

      if (vday[j]==0){
        Chat[,j] <- Chat[,j]-expbyglm
      }
    }
    MCE <- rep(0,length(vday))

    for (j in 1:length(vday)){
      MCE[j] <- sum(  (Chat[(timewt[,j] !=0),j]-tmp1[(timewt[,j] !=0),j])^2)
     }	
    print(cbind(vday,MCE))
    return(sum(MCE))
   }
  
# equn.15.min <- optim(start,newequn.15,method="Nelder-Mead",
#			control=list(maxit=200,trace=2),
#			Chat=Chat,lambdadx=lambdadx,w=w,v=v,vday=vday,
#			sigma2=sigma2,phi=phi,dv=dv,m=m,n=n,timewt=timewt, 
#			expbyglm=expbyglm)
 equn.15.min <- optimize(newequn.15,c(0.0000001,60),tol=0.0001,maximum=FALSE,
			Chat=Chat,lambdadx=lambdadx,w=w,v=v,vday=vday,
			sigma2=sigma2,phi=phi,dv=dv,m=128,n=128,timewt=timewt, 
			expbyglm=expbyglm,RATIO=RATIO)
  return(list(beta=equn.15.min$par[1]))
}




#theoritical covariance temporal function
cov.v.beta <- function(lambdadx,w,v,vday,sigma2,beta,phi,m,n){

#  res <-as.double(rep(0,length(vday)))

  .C("covequn",as.double(lambdadx),as.integer(128),as.integer(128),
   as.double(w),as.double(v),as.double(vday),as.integer(length(vday)),
     as.double(sigma2),as.double(beta),as.double(phi),res=double(length(vday)) )$res
}



#########covariate

exponential.cov <- function(x,beta)
{
  return(exp(-x/beta))
}

##################


find.fourier<-function(a,M,N,l1,l2,covfct) {
  # Finds the Fourier transform of the covariance function rho with
  # parameter a on the grid MxN, where the side lengths are supposed to be
  # l1 (x-direction) and l2 (y-direction).
  N2<-2*N
  M2<-2*M
  dd<-matrix(0,M2,N2)
  cc<-matrix(0,M2,N2)
  S<-matrix(0,M2,N2)
  delta1<-(l2/M)^2
  delta2<-(l1/N)^2
  for (i in 0:(M-1))
    for (j in 0:(N-1))
      dd[i+1,j+1]<-sqrt(delta1*i*i+delta2*j*j)
  for (i in M:(M2-1)) {
    for (j in 0:(N-1))
      dd[i+1,j+1]<-sqrt(delta1*(M2-i)*(M2-i)+delta2*j*j)
  }	
  for (i in 0:(M-1))
    for (j in N:(N2-1))
      dd[i+1,j+1]<-sqrt(delta1*i*i+delta2*(N2-j)*(N2-j))
  for (i in M:(M2-1))
    for (j in N:(N2-1))
      dd[i+1,j+1]<-sqrt(delta1*(M2-i)*(M2-i)+delta2*(N2-j)*(N2-j))
  covfct <- get(covfct,mode="function")
  cc <- covfct(dd,a)
  S<-fft(cc,inverse=F)
  return(S)
}


gauss.sim <- function(M,N,S) {
  # Calculates a gaussian random field, S is the Fourier transform of
  # the covariance matrix.
  N2<-2*N
  M2<-2*M
  X<-matrix(rnorm(N2*M2,0,1),M2,N2)
  X<-fft(X,inverse=F)
  A<-X*sqrt(S)
  W<-Re(fft(A,inverse=T))/(M2*N2)	
  G <- W[1:M,1:N]
  return(G)
}




sim.pois <- function(intensity,l1,l2) {
#  lambda <- sum(intensity)*l1*l2/prod(dim(intensity))
  lambda <- sum(intensity)
  #print(lambda)
  npts <- rpois(1,lambda)
  #print(paste("Simulating Poisson process with ",npts," points",sep=""))
  if (npts!=0) {
    spp.data <- matrix(0,ncol=2,nrow=npts)
    storage.mode(npts) <- "integer"
    N <- as.integer(length(intensity[,1]))
    M <- as.integer(length(intensity[1,]))
    intensity <- as.double(t(intensity))
    storage.mode(l1) <- "double"
    storage.mode(l2) <- "double"
    res <- matrix(.C("sim_pois",intensity,npts,l1,l2,N,M,spp=spp.data)$spp,
		  ncol=2,byrow=T)
  }
  else {
    res <- c(999)
  }
  return(res)
}


write.matrix <- function(x,file="",sep=" ",...)
{
  x <- as.matrix(x)
  p <- ncol(x)
  cat(format(t(x)),file=file,sep=c(rep(sep,p-1),"\n"),...)
}


rho2 <- function(x,beta,alpha,covfct)
{
  covfct <- get(covfct,mode="function")
  return(covfct(x,beta)^(2*alpha))
}

Aintegral <- function(epsilon,a0,alpha,beta,covfct)
{
  #return(integrate(rho2,epsilon,a0,100,500,0.01,beta=beta,alpha=alpha,
  return(integrate(rho2,epsilon,a0,100,500,0.0001,beta=beta,alpha=alpha,
     covfct=covfct)$value)
}

Bintegral <- function(g.dat,epsilon,a0,alpha,beta,covfct)
{
  covfct <- get(covfct,mode="function")
  g.dat <- g.dat[(g.dat[,2]>1),]
  c.fct <- log(g.dat[,2])
  r.fct <- covfct(g.dat[,1],beta)
  tmp <- c(epsilon,g.dat[1:(length(g.dat[,1])-1),1])
  dt <- g.dat[,1]-tmp
  return(sum((c.fct*r.fct)^alpha*dt))
}

minifct <- function(beta,g.dat,epsilon,a0,alpha,covfct)
{
  epsilonx <- epsilon[1]
  a0x <- a0[1]
  alphax <- alpha[1]
  covfctx <- as.character(covfct[1])
  return(Bintegral(g.dat,epsilonx,a0x,alphax,beta,covfctx)^2/
         Aintegral(epsilonx,a0x,alphax,beta,covfctx))
}

#lgcp.est <- function(g.dat,intensity,covfct,epsilon,a0,alpha)
lgcp.est <- function(g.dat,covfct,epsilon,a0,alpha)
{
  i1 <- g.dat[,1]>=epsilon
  i2 <- g.dat[,1]<=a0
  g.dat <- g.dat[i1&i2,]
  bhat <- optimise(minifct,c(0.01,1e6),0.01,1e6,maximum=TRUE,
               tol=.Machine$double.eps^0.25,g.dat,epsilon,
               a0,alpha,covfct)$objective
  shat <- (Bintegral(g.dat,epsilon,a0,alpha,bhat,covfct)/
         Aintegral(epsilon,a0,alpha,bhat,covfct))^(1/alpha)
#  mhat <- log(intensity)-log(sum(lambda))-shat/2
  return(list(beta=bhat,sigma=shat))
}

SIGMAFN <- function(g.dat,epsilon,a0,alpha,bhat,covfct){
shat <- (Bintegral(g.dat,epsilon,a0,alpha,bhat,covfct)/
         Aintegral(epsilon,a0,alpha,bhat,covfct))^(1/alpha)
return(shat)
}
