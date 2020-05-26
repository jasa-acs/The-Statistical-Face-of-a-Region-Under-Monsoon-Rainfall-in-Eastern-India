# This code is licensed under the MIT License

# Input: Starhull of complete rainfall regions over the years 2002-2012
# Output: Figure 4 of the paper

# Execution time: few seconds (Intel(R) Core(TM) i7-7700HQ CPU @ 2.80 GHz 2.81 GHz with RAM 16 GB)

# set directory
# setwd("/home")

# starhull contours over the years 2002-2012 
input_file="data/starhull_complete_contours"
list_file<- list.files(input_file, pattern = "\\.Rda$",full.names = TRUE) 

# read input
data1=NULL
for (i in 1:11){
dt3=NULL
dt3=readRDS(file=list_file[i]) # all contours of a year 
data1=cbind(data1,dt3)         # all complete larger contour function
}

# Normal transformation
s=1000                         
N11=ncol(data1)                # number of complete contours
r_u=(rank(data1))/(s*N11)      # ECDF
r_u=replace(r_u,which(r_u==1),0.9999)   
data2=matrix(qnorm(r_u),nrow=s) # Normal (random variable) transformation

# Skewness computation
s1=array(0,dim=1000)
s2=array(0,dim=1000)
for (i in 1:1000)
{
s1[i]=(mean((data1[i,]-mean(data1[i,]))^3))/((mean((data1[i,]-mean(data1[i,]))^2))^1.5)
s2[i]=(mean((data2[i,]-mean(data2[i,]))^3))/((mean((data2[i,]-mean(data2[i,]))^2))^1.5)
}

par(mfrow=c(1,1))
theta1=seq(0,2*pi,length.out = 1000)
plot(theta1,s1,col=1,lty=4,lwd=1.5,type="l",xlim=range(theta1),
     xaxt="n",ylim=range(-.7,3.5),cex.lab=0.8,ylab="Skewness",
     xlab="theta",cex.axis=.8)
axis(side=1, at=seq(0,2*pi,length.out=5),cex.axis=.8,
     labels=c(expression(paste(0,~(East))),expression(paste(pi/2, ~(North))),
              expression(paste(pi, ~(West))),expression(paste(3*pi/2, ~(South))),
              expression(paste(2*pi, ~(East)))))
lines(theta1,s2,col=4,lty=1,lwd=1.5)
legend(0.5,-0.08,c("Before transformation","After transformation"),col=c(1,4),cex=0.8,
       lty=c(4,1), lwd=c(1.5,1.5),box.lwd = 0,box.col = "white",bg = "white",horiz = TRUE)

