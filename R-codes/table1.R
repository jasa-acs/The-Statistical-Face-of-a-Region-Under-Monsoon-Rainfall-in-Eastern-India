# This code is licensed under the MIT License

# Input: Starhull approximated contours over the years 2002-2012
# Output: table 1 of the paper
# execution time: few seconds (Intel(R) Core(TM) i7-7700HQ CPU @ 2.80 GHz 2.81 GHz with RAM 16 GB)


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
table1=cumsum(pc_1)       # Cumulative percentage of variation  
print(table1)
