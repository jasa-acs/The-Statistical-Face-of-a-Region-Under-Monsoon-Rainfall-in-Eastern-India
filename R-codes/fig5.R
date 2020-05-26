# This code is licensed under the MIT License

# Input: Starhull approximated contours over the years 2002-2012
# Output: Figure 5 of the paper

# Execution time: 24 seconds (Intel(R) Core(TM) i7-7700HQ CPU @ 2.80 GHz 2.81 GHz with RAM 16 GB)


# set directory
# setwd("/home")

# starhull contours over the years 2002-2012 
input_file="data/starhull_complete_contours"
list_file<- list.files(input_file, pattern = "\\.Rda$",full.names = TRUE) 

# read input
data1=NULL
n1=NULL
for (i in 1:11){
  dt3=NULL
  dt3=readRDS(file=list_file[i]) # contours of i-th year 
  n1[i]=ncol(dt3)
  data1=cbind(data1,dt3)        
}

# Normal transformation
n11=cumsum(n1)
s=nrow(data1)
N11=ncol(data1)                   # Number of complete contours
r_u=(rank(data1))/(s*N11)         # Emperical CDF
r_u=replace(r_u,which(r_u==1),0.9999) 
data2=matrix(qnorm(r_u),nrow=s)   
data=data2
mean_x=rowMeans(data)             # Componentwise mean
cov1=array(0,dim=c(s,s))
cov1=(tcrossprod((data-mean_x)*(matrix(sqrt(N11), nrow = s, 
                ncol = ncol(data), byrow = TRUE))))     

# Eigen analysis of complete data 
h1=eigen(cov1)$values  
hv=eigen(cov1)$vectors  
pc_1=((h1[1:12])/sum(h1))*100;pc_1 # Percentage of the variation explained by 12 PCS

# componentwise plot in linear scale (re-transformed scale) 
par(mfrow=c(1,1),mar = c(.1,.1,.5,.5),oma = c(4,4,.1,.1), xpd = NA)            
theta=seq(0,2*pi,length.out=(s+1))[-(s+1)]
plot(0,0,type="l",xaxt='n',xlim=c(0,2*pi),ylim=c(-0.055,0.05),main="",
     cex.axis=0.8,cex.lab=0.8,cex=1,xlab="Angle",ylab="Value of PC function")  
axis(side=1, at=seq(0,2*pi,length.out=5),cex.axis=.8,
     labels=c(expression(paste(0,~(East))),expression(paste(pi/2, ~(North))),
              expression(paste(pi, ~(West))),expression(paste(3*pi/2, ~(South))),
              expression(paste(2*pi, ~(East)))))

# plot of first 10 eigen functions
M=10
for (j in 1:M)
{
  r2=hv[,j]
 if (j==1){
   lines(theta,-r2,type="l",ylim=c(-0.06,0.06),cex=0.9,lty=j,col=j,cex.main=1,lwd=1,
         xlab="Angle",ylab="Value of PC function")  
 }else {
   lines(theta,r2,type="l",ylim=c(-0.06,0.06),cex=0.9,lty=j,col=j,lwd=1,
         cex.main=1,xlab="Angle",ylab="Value of PC function")  
}
   }
legend(-0.15,-0.05,c("PC1", "PC2", "PC3","PC4", "PC5","PC6", "PC7", "PC8","PC9","PC10"), lwd=1, 
       cex=0.6, border=NA,horiz = TRUE,lty=c(1:10),col=c(1:10))


