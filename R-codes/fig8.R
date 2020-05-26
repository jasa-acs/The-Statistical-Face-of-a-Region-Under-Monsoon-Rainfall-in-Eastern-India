# This code is licensed under the MIT License

# Input: Starhull over the years 2002-2012
# Output: Figure 8 of the paper


# load library 
library(sp)
library(geosphere)
library(survival)

# set directory
# setwd("/home")

# input data
input_file1 <- "data/starhull_complete_contours"
input_file2 <- "data/area_yearwise"
input_file3 <- "data/censor_status_yearwise"
file1 <- list.files(input_file1, pattern = "\\.Rda$",full.names = TRUE) 
file2 <- list.files(input_file2, pattern = "\\.Rda$",full.names = TRUE)   
file3 <- list.files(input_file3, pattern = "\\.Rda$",full.names = TRUE)   
data1 <- NULL
n1 <- NULL
for (j in 1:11){
  ar <- NULL
  ar <- readRDS(file=file1[j])
  n1[j] <- ncol(ar)
  data1 <- cbind(data1,ar) 
}
data_total <- data1   
n11 <- c(0,cumsum(n1))
s <- 1000                  
N11 <- ncol(data1)   
r_u <- (rank(data1))/(s*N11)               # Emperical DF
r_u <- replace(r_u,which(r_u==1),0.9999)  
data_tarns <- matrix(qnorm(r_u),nrow=s)   # Normal (random variable) transformation

# eigen function sign assignment for plot
sign1 <- array(0,dim=c(9,11))
x1 <- c(1,1,1,1,1,1,1,1,1,1,1)
x2 <- c(-1,1,1,1,-1,-1,1,-1,1,1,1)
x3 <- -c(1,1,-1,-1,-1,-1,1,-1,-1,-1,-1)
x4 <- c(1,-1,1,-1,1,1,1,-1,-1,1,-1)
x5 <- c(-1,1,-1,-1,1,-1,-1,-1,-1,-1,-1)
x6 <- c(1,1,-1,1,1,1,1,1,-1,1,1)
x7 <- c(-1,1,1,-1,-1,-1,1,1,1,1,1)
x8 <- c(1,-1,-1,1,-1,-1,1,-1,-1,-1,-1)
x9 <- c(-1,-1,-1,1,-1,-1,-1,-1,-1,-1,-1)
sign1 <- rbind(x1,x2,x3,x4,x5,x6,x7,x8,x9)

# plot configuration
m <- matrix(c(1,2,3,4,5,6,7,8,9,10,10,10),nrow = 4,ncol = 3,byrow = TRUE)
layout(mat = m,heights = c(2,2,2))
par(mar = c(0.1, 0.1, 0.1, 0.1), oma = c(4,4,0,0), xpd = NA)         

for (k in 1:9){                               # plot one panel (one year) 
  plot(0,0, yaxt='n',xaxt='n', axes=FALSE, frame.plot=TRUE,ann=FALSE,col="white",type="l",ylim=c(-0.08,0.08),xlim=c(0,2*pi),cex.main=1)  
  y=seq(-0.04,0.04,length.out=3)
  if (k>=7) {axis(side=1, at=seq(0,2*pi,length.out=5),cex.axis=1,
                  labels=c(expression(0),expression(pi/2),expression(pi),expression(3*pi/2),expression(2*pi)))}
  if (k==1||k==4||k==7){axis(side = 2,cex.axis=1,at=y, labels = if (k==1||k==4||k==7) levels(y) else FALSE)}

  cov1 <- array(0,dim=c(s,s))
  mean_x <- array(0,dim=s)
  for (m in 1:11)                                # one curve (PC function) in each panel
  {
    ar <- readRDS(file=file2[m])
    st <- readRDS(file=file3[m])
    dt4 <- intersect(which(ar>=200),which(ar<=13500)) # All contours in the bound
    
    dt1 <- which(st==1)                    # censored regions
    dt2 <- which(st==2)                    # complete regions
    
    dt11 <- intersect(dt1,dt4)             # (> 200 sq km && <13500) censored 
    dt22 <- intersect(dt2,dt4)             # (> 200 sq km && <13500) complete 
    N1 <- length(dt22)                      
    N2 <- length(dt11)                       
    
    time <- NULL
    censor <- NULL
    time <- ar[dt4]                        # contour lies in the given threshold
    censor <- st[dt4]
    censor[censor < 1.5] <- 0              # assign censored  0
    censor[censor > 1.5] <- 1              # assign complete  1
    
    time_ind <- order(time)
    area_order <- time[time_ind]
    censor_order <- censor[time_ind]
    area_status <- cbind(area_order,censor_order)
    
    c_max <- max(which(area_status[,2]==1))
    area_status1 <- area_status[1:c_max,]
    time <- area_status1[,1]
    censor <- area_status1[,2]
    
    mini.surv <- survfit(Surv(time, censor)~ 1)    # Kaplan-Meier computation
    T1 <- mini.surv$time                           # increasing time period 
    F1 <- mini.surv$surv                           # survial function
    F2 <- c(1,F1[-length(F1)])-F1
    
    F3 <- cbind(T1,F2)                             # (index,prob. mass)
    F4 <- array(0,dim = N1)
    F4 <- c(F3[which(F3[,2]>0),2])                 # Weights
    
    order_com <- order(ar[dt22])
    x0 <- order_com                                # Order the complete contours
    data_star <- data_tarns[,(n11[m]+1):n11[m+1]]  # data corresponding to one year
    N1 <- length(x0)
    mean_x <- c(data_star[,x0]%*%F4)               # Weighted mean
    cov1 <- (tcrossprod((data_star[,x0]-mean_x)*(matrix(sqrt(F4), nrow = s, ncol = N1, byrow = TRUE)))) 
    h1 <- eigen(cov1)$values                       
    hv <- eigen(cov1)$vectors                      
    
    theta <- seq(0,2*pi,length.out= (s+1))[-(s+1)]
    pc_1 <- ((h1[1:20])/sum(h1))*100;pc_1 
    r2 <- hv[,k]*sign1[k,m]                       # Use sign matrix
    lines(theta,r2,type="l",col=(m+2),ylim=c(-0.08,0.08),
          main="PC",cex.main=0.8,lty=(m+2))
  }
}

plot(1, type = "n", axes=FALSE, xlab="", ylab="")
plot_colors <- seq(1,11,1)+2
legend(x = "top",inset =.7,
       legend = c( "02", "03","04", "05","06", "07", "08","09", "10","11","12"), 
       col=plot_colors, lwd=1, cex=0.9, horiz = TRUE,lty=(c(1:11)+2))
title(xlab = "Direction (from due East)",
      ylab = "Value of PC function",
      outer = TRUE, line = 2)

