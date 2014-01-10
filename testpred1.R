source(paste(UTILITYDIR,"parestimate3.r",sep=""))
M <- 128
N <- 128


#when use read.table, read in as a list
INTENSITY3 <- read.table(paste(DATADIR,"INTENSITY3.txt",sep=""),header=F)
isinw <- read.table(paste(DATADIR,"isinw.txt",sep=""),header=F)

#Create sommthdens.in from INTENSITY3.txt.
#Prepare the use of multimala.in
abc <- matrix(0,2*M,2*N)
abc[1:M,1:N] <- as.matrix(INTENSITY3)
write.matrix(abc,paste(UTILITYDIR,"smoothdens.in",sep=''))


#Create isinw.in from isinw.txt
#Prepare for multimala.in
abc <- matrix(0,2*M,2*N)
abc[1:M,1:N] <- as.matrix(isinw)
write.matrix(abc,paste(UTILITYDIR,"isinw.in",sep=''))

#Phi is the spatio parameter.
#Doing the folloing step generates S.in
#Current S.in is based on Phi<-6.525996
S_find.fourier(Phi,M,N,l1,l2,"exponential.cov")
S1_t(Re(S))
if (length(S1[S1<0])>0)
print("Warning: Negative definite covariance for extended field.")
write.matrix(S1,paste(UTILITYDIR,"S.in",sep=''))

 #CREATE DATA-3 cols:X/Y/ONSETDATE.  X/Y in Km.

temp2 <- sppdata(DATA,c(402,491),c(87,176))
write(t(temp2),file="~/AEGISS/LGCP/UTILITY/spp1",ncolumns=3)
