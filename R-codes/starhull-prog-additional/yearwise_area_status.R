# This code is licensed under the MIT License

# Input: contours and missing status of contours 
# over the years 2002-2012.
# Output: area of rainfall region and its missing status

# load libraries
library(sp)
library(geosphere)


# function to compute the area under the rainfall regions
func3 <-  function(xx,yy,zz){
  p1 <-  as.data.frame(xx)
  p2 <-  as.data.frame(yy)
  p11 <-  as.data.frame(zz)
  p3 <-  p1[-which(p1[,1]==0),]           
  p5 <-  p2[-which(p2[,1]==0),]
  p111 <-  p11[-which(p11[,1]==0),]
  p4 <-  cumsum(p5)                      
  N <-  length(p4)
  p4 <-  c(0,p4)
  contour_prj <-  list()
  area1 <-  array(0,dim =  N)
  q <-  list()
  ll <-  NULL
  for (i in 1:N)
  {
    obj_1 <-  NULL
    obj_1 <-  p3[(p4[i]+1):p4[i+1],]      # Shape objects_1
    if (nrow(obj_1)>= 4){    
      m00 <-  NULL
      m00 <- data.frame(x =  obj_1[,1], y =  obj_1[,2])
      area1[i] <-  areaPolygon(m00)/(10^6)   # Compute the area of a polygon incoordinates (longitude/latitude)
      }
  }
  return(cbind(area1,p111))
}

# In each the following block data corresponding to each year
# taken as input of the above function

# year 2002
h_6 <-  func3(read.table("data/rf_data/2002/M67_cluster_points"),
       read.table("data/rf_data/2002/M67_cluster_lengths"),
       read.table("data/rf_data/2002/M67_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2002/M89_cluster_points"),
       read.table("data/rf_data/2002/M89_cluster_lengths"),
       read.table("data/rf_data/2002/M89_cen_spell"))
h2 <-  rbind(h_6,h_7)
ar1 <-  c(h_6[,1],h_7[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2])   # Status fo the region (censored or complete)

# write results
saveRDS(ar1,file <-  "data/area_yearwise/2002.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2002.Rda")

# year 2003
h_6 <-  func3(read.table("data/rf_data/2003/M67_cluster_points"),
       read.table("data/rf_data/2003/M67_cluster_lengths"),
       read.table("data/rf_data/2003/M67_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2003/M89_cluster_points"),
       read.table("data/rf_data/2003/M89_cluster_lengths"),
       read.table("data/rf_data/2003/M89_cen_spell"))
h3 <-  rbind(h_6,h_7)
ar1 <-  NULL
st <-  NULL
ar1 <-  c(h_6[,1],h_7[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2])   # Status fo the region (censored or complete)
saveRDS(ar1,file <-  "data/area_yearwise/2003.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2003.Rda")

# year 2004
h_6 <-  func3(read.table("data/rf_data/2004/M6_cluster_points"),
       read.table("data/rf_data/2004/M6_cluster_lengths"),
       read.table("data/rf_data/2004/M6_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2004/M7_cluster_points"),
       read.table("data/rf_data/2004/M7_cluster_lengths"),
       read.table("data/rf_data/2004/M7_cen_spell"))
h_8 <-  func3(read.table("data/rf_data/2004/M8_cluster_points"),
       read.table("data/rf_data/2004/M8_cluster_lengths"),
       read.table("data/rf_data/2004/M8_cen_spell"))
h_9 <-  func3(read.table("data/rf_data/2004/M9_cluster_points"),
       read.table("data/rf_data/2004/M9_cluster_lengths"),
       read.table("data/rf_data/2004/M9_cen_spell"))
h4 <-  rbind(h_6,h_7,h_8,h_9)
ar1 <-  NULL
st <-  NULL
ar1 <-  c(h_6[,1],h_7[,1],h_8[,1],h_9[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2],h_8[,2],h_9[,2])   # Status fo the region (censored or complete)
saveRDS(ar1,file <-  "data/area_yearwise/2004.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2004.Rda")

# year 2005
h_6 <-  func3(read.table("data/rf_data/2005/M6_cluster_points"),
       read.table("data/rf_data/2005/M6_cluster_lengths"),
       read.table("data/rf_data/2005/M6_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2005/M7_cluster_points"),
       read.table("data/rf_data/2005/M7_cluster_lengths"),
       read.table("data/rf_data/2005/M7_cen_spell"))
h_8 <-  func3(read.table("data/rf_data/2005/M8_cluster_points"),
       read.table("data/rf_data/2005/M8_cluster_lengths"),
       read.table("data/rf_data/2005/M8_cen_spell"))
h_9 <-  func3(read.table("data/rf_data/2005/M9_cluster_points"),
       read.table("data/rf_data/2005/M9_cluster_lengths"),
       read.table("data/rf_data/2005/M9_cen_spell"))
h5 <-  rbind(h_6,h_7,h_8,h_9)
ar1 <-  NULL
st <-  NULL
ar1 <-  c(h_6[,1],h_7[,1],h_8[,1],h_9[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2],h_8[,2],h_9[,2])   # Status fo the region (censored or complete)
saveRDS(ar1,file <-  "data/area_yearwise/2005.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2005.Rda")

# year 2006
h_6 <-  func3(read.table("data/rf_data/2006/M6_cluster_points"),
       read.table("data/rf_data/2006/M6_cluster_lengths"),
       read.table("data/rf_data/2006/M6_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2006/M7_cluster_points"),
       read.table("data/rf_data/2006/M7_cluster_lengths"),
       read.table("data/rf_data/2006/M7_cen_spell"))
h_8 <-  func3(read.table("data/rf_data/2006/M8_cluster_points"),
       read.table("data/rf_data/2006/M8_cluster_lengths"),
       read.table("data/rf_data/2006/M8_cen_spell"))
h_9 <-  func3(read.table("data/rf_data/2006/M9_cluster_points"),
       read.table("data/rf_data/2006/M9_cluster_lengths"),
       read.table("data/rf_data/2006/M9_cen_spell"))
h6 <-  rbind(h_6,h_7,h_8,h_9)
ar1 <-  NULL
st <-  NULL
ar1 <-  c(h_6[,1],h_7[,1],h_8[,1],h_9[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2],h_8[,2],h_9[,2])   # Status fo the region (censored or complete)
saveRDS(ar1,file <-  "data/area_yearwise/2006.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2006.Rda")

# year 2007
h_6 <-  func3(read.table("data/rf_data/2007/M6_cluster_points"),
       read.table("data/rf_data/2007/M6_cluster_lengths"),
       read.table("data/rf_data/2007/M6_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2007/M7_cluster_points"),
       read.table("data/rf_data/2007/M7_cluster_lengths"),
       read.table("data/rf_data/2007/M7_cen_spell"))
h_8 <-  func3(read.table("data/rf_data/2007/M8_cluster_points"),
       read.table("data/rf_data/2007/M8_cluster_lengths"),
       read.table("data/rf_data/2007/M8_cen_spell"))
h_9 <-  func3(read.table("data/rf_data/2007/M9_cluster_points"),
       read.table("data/rf_data/2007/M9_cluster_lengths"),
       read.table("data/rf_data/2007/M9_cen_spell"))
h7 <-  rbind(h_6,h_7,h_8,h_9)
ar1 <-  NULL
st <-  NULL
ar1 <-  c(h_6[,1],h_7[,1],h_8[,1],h_9[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2],h_8[,2],h_9[,2])   # Status fo the region (censored or complete)
saveRDS(ar1,file <-  "data/area_yearwise/2007.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2007.Rda")

# year 2008
h_6 <-  func3(read.table("data/rf_data/2008/M6_cluster_points"),
       read.table("data/rf_data/2008/M6_cluster_lengths"),
       read.table("data/rf_data/2008/M6_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2008/M7_cluster_points"),
       read.table("data/rf_data/2008/M7_cluster_lengths"),
       read.table("data/rf_data/2008/M7_cen_spell"))
h_8 <-  func3(read.table("data/rf_data/2008/M8_cluster_points"),
       read.table("data/rf_data/2008/M8_cluster_lengths"),
       read.table("data/rf_data/2008/M8_cen_spell"))
h_9 <-  func3(read.table("data/rf_data/2008/M9_cluster_points"),
       read.table("data/rf_data/2008/M9_cluster_lengths"),
       read.table("data/rf_data/2008/M9_cen_spell"))
h8 <-  rbind(h_6,h_7,h_8,h_9)
ar1 <-  NULL
st <-  NULL
ar1 <-  c(h_6[,1],h_7[,1],h_8[,1],h_9[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2],h_8[,2],h_9[,2])   # Status fo the region (censored or complete)
saveRDS(ar1,file <-  "data/area_yearwise/2008.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2008.Rda")

# year 2009
h_6 <-  func3(read.table("data/rf_data/2009/M6_cluster_points"),
       read.table("data/rf_data/2009/M6_cluster_lengths"),
       read.table("data/rf_data/2009/M6_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2009/M7_cluster_points"),
       read.table("data/rf_data/2009/M7_cluster_lengths"),
       read.table("data/rf_data/2009/M7_cen_spell"))
h_8 <-  func3(read.table("data/rf_data/2009/M8_cluster_points"),
       read.table("data/rf_data/2009/M8_cluster_lengths"),
       read.table("data/rf_data/2009/M8_cen_spell"))
h_9 <-  func3(read.table("data/rf_data/2009/M9_cluster_points"),
       read.table("data/rf_data/2009/M9_cluster_lengths"),
       read.table("data/rf_data/2009/M9_cen_spell"))
h9 <-  rbind(h_6,h_7,h_8,h_9)
ar1 <-  NULL
st <-  NULL
ar1 <-  c(h_6[,1],h_7[,1],h_8[,1],h_9[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2],h_8[,2],h_9[,2])   # Status fo the region (censored or complete)
saveRDS(ar1,file <-  "data/area_yearwise/2009.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2009.Rda")

# year 2010
h_6 <-  func3(read.table("data/rf_data/2010/M6_cluster_points"),
       read.table("data/rf_data/2010/M6_cluster_lengths"),
       read.table("data/rf_data/2010/M6_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2010/M7_cluster_points"),
       read.table("data/rf_data/2010/M7_cluster_lengths"),
       read.table("data/rf_data/2010/M7_cen_spell"))
h_8 <-  func3(read.table("data/rf_data/2010/M8_cluster_points"),
       read.table("data/rf_data/2010/M8_cluster_lengths"),
       read.table("data/rf_data/2010/M8_cen_spell"))
h_9 <-  func3(read.table("data/rf_data/2010/M9_cluster_points"),
       read.table("data/rf_data/2010/M9_cluster_lengths"),
       read.table("data/rf_data/2010/M9_cen_spell"))
h10 <-  rbind(h_6,h_7,h_8,h_9)
ar1 <-  NULL
st <-  NULL
ar1 <-  c(h_6[,1],h_7[,1],h_8[,1],h_9[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2],h_8[,2],h_9[,2])   # Status fo the region (censored or complete)
saveRDS(ar1,file <-  "data/area_yearwise/2010.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2010.Rda")

# year 2011
h_6 <-  func3(read.table("data/rf_data/2011/M6_cluster_points"),
       read.table("data/rf_data/2011/M6_cluster_lengths"),
       read.table("data/rf_data/2011/M6_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2011/M7_cluster_points"),
       read.table("data/rf_data/2011/M7_cluster_lengths"),
       read.table("data/rf_data/2011/M7_cen_spell"))
h_8 <-  func3(read.table("data/rf_data/2011/M8_cluster_points"),
       read.table("data/rf_data/2011/M8_cluster_lengths"),
       read.table("data/rf_data/2011/M8_cen_spell"))
h_9 <-  func3(read.table("data/rf_data/2011/M9_cluster_points"),
       read.table("data/rf_data/2011/M9_cluster_lengths"),
       read.table("data/rf_data/2011/M9_cen_spell"))
h11 <-  rbind(h_6,h_7,h_8,h_9)
ar1 <-  NULL
st <-  NULL
ar1 <-  c(h_6[,1],h_7[,1],h_8[,1],h_9[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2],h_8[,2],h_9[,2])   # Status fo the region (censored or complete)
saveRDS(ar1,file <-  "data/area_yearwise/2011.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2011.Rda")

# year 2012
h_6 <-  func3(read.table("data/rf_data/2012/M6_cluster_points"),
       read.table("data/rf_data/2012/M6_cluster_lengths"),
       read.table("data/rf_data/2012/M6_cen_spell"))
h_7 <-  func3(read.table("data/rf_data/2012/M7_cluster_points"),
       read.table("data/rf_data/2012/M7_cluster_lengths"),
       read.table("data/rf_data/2012/M7_cen_spell"))
h_8 <-  func3(read.table("data/rf_data/2012/M8_cluster_points"),
       read.table("data/rf_data/2012/M8_cluster_lengths"),
       read.table("data/rf_data/2012/M8_cen_spell"))
h_9 <-  func3(read.table("data/rf_data/2012/M9_cluster_points"),
       read.table("data/rf_data/2012/M9_cluster_lengths"),
       read.table("data/rf_data/2012/M9_cen_spell"))
h12 <-  rbind(h_6,h_7,h_8,h_9)
ar1 <-  NULL
st <-  NULL
ar1 <-  c(h_6[,1],h_7[,1],h_8[,1],h_9[,1])  # Area of regions
st <-  c(h_6[,2],h_7[,2],h_8[,2],h_9[,2])   # Status fo the region (censored or complete)
saveRDS(ar1,file <-  "data/area_yearwise/2012.Rda")
saveRDS(st,file <-  "data/censor_status_yearwise/2012.Rda")

ar <-  c(h2[,1],h3[,1],h4[,1],h5[,1],h6[,1],h7[,1],h8[,1],h9[,1],h10[,1],h11[,1],h12[,1])  # Area of regions
st <-  c(h2[,2],h3[,2],h4[,2],h5[,2],h6[,2],h7[,2],h8[,2],h9[,2],h10[,2],h11[,2],h12[,2])   # Status fo the region (censored or complete)

# write combined result
saveRDS(ar,file <-  "data/area_status_combined/area.Rda")
saveRDS(st,file <-  "data/area_status_combined/status.Rda")
