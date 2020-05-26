# This code is licensed under the MIT License

# Input: Starhull approximated contours of rainfall regions for the years 2002-2012
# Output: Figure 9 of the paper 


# load libraries
library(sp)
library(geosphere)
library(survival)

# set directory
# setwd("/home")

# input files
file1<- list.files("data/area_yearwise", pattern = "\\.Rda$",full.names = TRUE)   
file2<- list.files("data/censor_status_yearwise", pattern = "\\.Rda$",full.names = TRUE)   
file3<- list.files("data/starhull_complete_contours", pattern = "\\.Rda$",full.names = TRUE) 
data1=NULL
n1=NULL
for (j in 1:11){
  dt3=NULL
  dt3=readRDS(file=file3[j])       # all contours of a year 
  n1[j]=ncol(dt3)
  data1=cbind(data1,dt3)
}
data_total=data1   
n11=c(0,cumsum(n1))
s=1000
N11=ncol(data1)                       # Number of complete contours
r_u=(rank(data1))/(s*N11)             # ECDF
r_u=replace(r_u,which(r_u==1),0.9999) 
data_tarns=matrix(qnorm(r_u),nrow=s)   # Normal (random variable) transformation

# Yearly eigen decomposition
pc_1=array(0,dim=c(11,12))
plot(x=1,y=1,type='n',log="y",yaxt="n",ylim=range(1,1300),xlim=c(0.1,20),ylab="Eigenvalue",xlab="Component index",lwd=0,cex.lab=0.8,cex.axis=0.8,cex.main=0.8)                    # Plot of the eigen values
axis(2, at=c(2,5,20,100,1000),cex.lab=0.8,cex.axis=0.8,cex.main=0.8)
for (m in 1:11)
  {
    dt3=readRDS(file=file1[m])
    st=readRDS(file=file2[m])
    dt4=intersect(which(dt3>=200),which(dt3<=13500))
    
    dt1=which(st==1)                    # censored regions
    dt2=which(st==2)                    # complete regions
    dt11=intersect(dt1,dt4)             # area (> 200 sq km && <13500) censored spell
    dt22=intersect(dt2,dt4)             # area (> 200 sq km && <13500) complete spell (except one or two)
    time1=NULL
    censor1=NULL
    time1=dt3[dt4]   
    censor1=st[dt4]  
    censor1[censor1 < 1.5]=0            # Censored 0
    censor1[censor1 > 1.5]=1            # Complete 1
    
    time_ind=order(time1)
    area_order=time1[time_ind]
    censor_order=censor1[time_ind]
    area_status=cbind(area_order,censor_order)
    
    c_max=max(which(area_status[,2]==1))
    area_status1=area_status[1:c_max,]
    time=NULL
    censor=NULL
    time=rank(area_status1[,1])
    censor=area_status1[,2]
    mini.surv <- survfit(Surv(time, censor)~ 1) # Kaplan-Meier computation
    T1=mini.surv$time                           # Increasing time period 
    F1=mini.surv$surv                          # Estimate of survial function
    F2=c(1,F1[-length(F1)])-F1
    F3=cbind(T1,F2)                            # (index,prob. mass)
    F4=c(F3[which(F3[,2]>0),2])                # Weights
    
    order_com=order(dt3[dt22])
    x0=order_com
    N1=length(x0)
    data_star=data_tarns[,(n11[m]+1):n11[m+1]]  # Particular year (2011) transformed dat
    mean_x=c(data_star[,x0]%*%F4)               # Weighted mean
    cov1=(tcrossprod((data_star[,x0]-mean_x)*(matrix(sqrt(F4), nrow = s, ncol = N1, byrow = TRUE))))  #  This is not covarince
    h1=eigen(cov1)$values  
    pc_1[m,]=((h1[1:12])/sum(h1))*100       # Percentage of the variation explained by the PCS
    lines(h1[1:20],type = "l",log="y",col=m,lty=m,lwd=1,cex.lab=0.8,cex.axis=0.8,cex.main=0.8)                    # Plot of the eigen values
     }
print(pc_1)  # table 2
legend("bottomright", legend = c( "02", "03","04", "05","06", "07", "08","09", "10","11","12"), 
col=seq(1,11,1),lty=seq(1,11,1),lwd=1, cex=0.55, horiz = TRUE)

# plot(log(h1[1:20]),type = "l",ylab="Eigenvalue",xlab="PC number",col=m,lty=m,lwd=1,cex.lab=0.8,cex.axis=0.8,cex.main=0.8)                    # Plot of the eigen values

