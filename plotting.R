library(shapelib)
library(splancs)
###define directory
PROGRAMDIR <- "~/AEGISS/LGCP/PROGRAM/"
DATADIR <- "~/AEGISS/LGCP/DATATOUSE/"
UTILITYDIR <- "~/AEGISS/LGCP/UTILITY/"
SIMDATADIR <- "~/AEGISS/LGCP/PROGRAM/SIMDATA/"
SIM0319 <- "~/AEGISS/LGCP/PROGRAM/SIMDATA/SIM0319/"
SIMOUTPUT <- "~/AEGISS/LGCP/UTILITY/SIMOUTPUT/"
source(paste(UTILITYDIR,"parestimate3.r",sep=""))


isinw <-read.table(paste(DATADIR,"isinw.txt",sep=""),header=F)
isinw <- as.matrix(isinw)

INTENSITY3 <- read.table(paste(DATADIR,"INTENSITY3.txt",sep=""),header=F)
INTENSITY3 <- as.matrix(INTENSITY3)

xgrid_seq(402000+(89000/128/2),491000-(89000/128/2),l=128)/1000
ygrid_seq(87000+(89000/128/2),176000-(89000/128/2),l=128)/1000

par(mfrow=c(2,2))
exceedc2<-matrix(scan(paste(SIMOUTPUT,"exceed5_4_2.out",sep="")),ncol=2*128, byrow=T)
exceedc2<-exceedc2[1:128,1:128]
tmp2_exceedc2*isinw/10000
tmp2[isinw==0]_NA
temp2_tmp2
temp2[is.na(temp2)] <-4
temp2[(tmp2<0.5)]_3
temp2[(tmp2 >= 0.5) & (tmp2<0.9)]_2
temp2[(tmp2>=0.9) & (tmp2<0.99)]_1
temp2[(tmp2>=0.99)]_0

exceedc4<-matrix(scan(paste(SIMOUTPUT,"exceed5_4_4.out",sep="")),ncol=2*128, byrow=T)
exceedc4<-exceedc4[1:128,1:128]
tmp4_exceedc4*isinw/10000
tmp4[isinw==0]_NA
temp4_tmp4
temp4[is.na(temp4)] <-4
temp4[(tmp4<0.5)]_3
temp4[(tmp4 >= 0.5) & (tmp4<0.9)]_2
temp4[(tmp4>=0.9) & (tmp4<0.99)]_1
temp4[(tmp4>=0.99)]_0

exceedc8<-matrix(scan(paste(SIMOUTPUT,"exceed5_4_8.out",sep="")),ncol=2*128, byrow=T)
exceedc8<-exceedc8[1:128,1:128]
tmp8_exceedc8*isinw/10000
tmp8[isinw==0]_NA
temp8_tmp8
temp8[is.na(temp8)] <-4
temp8[(tmp8<0.5)]_3
temp8[(tmp8 >= 0.5) & (tmp8<0.9)]_2
temp8[(tmp8>=0.9) & (tmp8<0.99)]_1
temp8[(tmp8>=0.99)]_0
#grey 0 is black, 1 is white

par(mfrow=c(3,2))
par(pty="s")

 image(xgrid*1000,ygrid*1000,t(temp2),xlim=c(400000,492000),ylim=c(87000,177000),zlim=c(0,4),xlab='Easting',ylab='Northing',col=c("red","blue","yellow","limegreen","white"))

legend(405000,170000,c("NA","0<=p<0.5","0.5<=p<0.9","0.9<=p<0.99","0.99<=p"),col=c("white","limegreen","yellow","blue","red"),fill=c("white","limegreen","yellow","blue","red"))


plot(c(402000,491000),c(87000,176000),type='n',xlab='Easting',ylab='Northing')
shapemap('/home/sut3/AEGISS/BOUNDRY/HAMP/edb25_c.shp',add=T)






