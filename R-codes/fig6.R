# This code is licensed under the MIT License

# Input: Starhull, area and missing status of regions 
# over the years 2002-2012.
# Output: Figure 6 of the paper.
# Execution time: few seconds (Intel(R) Core(TM) i7-7700HQ CPU @ 2.80 GHz 2.81 GHz with RAM 16 GB)


# load libraries
library(sp)
library(geosphere)
library(survival)

# set directory
# setwd("/home")

# input files
input_file1="data/starhull_complete_contours"    # starhull for 2002-2012
input_file2="data/area_status_combined/area.Rda" # area of the rainfall regions
input_file3="data/area_status_combined/status.Rda" # complete (2) or censored (1)
  
# read input
list_file= list.files(input_file1, pattern = "\\.Rda$",full.names = TRUE) 
ar=readRDS(file=input_file2)
st= readRDS(file=input_file3)

data1=NULL
n1=NULL
for (j in 1:11){
  dt=NULL
  dt=readRDS(file=list_file[j]) # All contours of a year 
  n1[j]=ncol(dt)
  data1=cbind(data1,dt) # all complete larger contour function
}
data_total=data1   # no: 11140
n11=cumsum(n1)
s=1000        

# normalization
N11=ncol(data1)   # Number of complete contours
r_u=(rank(data1))/(s*N11) # ECDF
r_u=replace(r_u,which(r_u==1),0.9999)   ### Check this line again
data2=matrix(qnorm(r_u),nrow=s)   # Normal (random variable) transformation
data_star=data2

# Computation of Kaplan-Meier weight
dt4=intersect(which(ar>=200),which(ar<=13500)) # All contours in the bound
dt1=which(st==1)  # Indices of all censored regions
dt2=which(st==2)  # Indices of all complete regions

dt11=intersect(dt1,dt4)  # Indices of larger (> 200 sq km && <13500) censored spell
dt22=intersect(dt2,dt4)  # Indices of larger (> 200 sq km && <13500) complete spell 

N1=length(dt22)  # Number of complete observations
N2=length(dt11) # Number of censored observations

time1=NULL
censor1=NULL
time1=ar[dt4]    # Region (Area) which lies in the given threshold
censor1=st[dt4]
censor1[censor1 < 1.5]=0   #Censored is 0
censor1[censor1 > 1.5]=1  # Complete observations is 1

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

mini.surv <- survfit(Surv(time, censor)~ 1) # Kaplan-Meier fit
T1=mini.surv$time  # Increasing time period 
F1=mini.surv$surv  # Estimate of survial function

F2=c(1,F1[-length(F1)])-F1
F3=cbind(T1,F2)               # (index,prob. mass)
F4=c(F3[which(F3[,2]>0),2])   # Weights

order_com=order(ar[dt22])
x0=order_com # Order the complete contours 
N1=length(x0)

mean_x=c(data_star[,x0]%*%F4)   # Weighted mean
cov1=(tcrossprod((data_star[,x0]-mean_x)*(matrix(sqrt(F4), nrow = s, ncol = N1, byrow = TRUE))))  
h1=eigen(cov1)$values 
hv=eigen(cov1)$vectors
h2=eigen(cov1)$values[1:20]   # first 20 eigen values

par(mfrow=c(1,1),mar = c(.1,.1,.5,.5),oma = c(4,4,.5,.5), xpd = NA)            # allow content to protrude into outer margin (and beyond)
plot(h2,type = "b",ylab="Eigenvalue",xlab="Principal Component number",lwd=1,cex.lab=0.9,cex.axis=0.9,cex.main=0.9)                    # Plot of the eigen values

pc_1=((h1[1:20])/sum(h1))*100;pc_1 # Percentage of the variation explained by the PCS
print(cumsum(pc_1))

