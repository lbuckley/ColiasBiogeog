
#========================================
#LOAD LIBRARIES

memory.size(max=TRUE)
memory.limit(size = 4095)

library(zoo)
library(chron) #convert dates
library(RAtmosphere)
library(msm)
library(MASS)
library(SDMTools)
library(gdata) #for converting to matrix
library(truncnorm)
library(ncdf)
library(chron)
library(foreach)
library(abind) #combine matrices into array
library(doParallel)
registerDoParallel(cl=4)

library(reshape2)
library(ggplot2)
library(grid)
library(data.table)
library(raster)
library(ks)
library(ncdf4)

#source functions
source("ColiasFunctions.R" )

#========================================

#CHOOSE PROJECTION
proj.k=2 #1: bcc-csm1-1.1.rcp60, 2: ccsm4.1.rcp60, 3: gfdl-cm3.1.rcp60
projs=c("bcc-csm","ccsm4","gfdl")

#==================================
#LOAD DATA

#Read DCP data
# http://gdo-dcp.ucllnl.org/downscaled_cmip_projections/dcpInterface.html

setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/Butterflies/Plasticity/ClimateData/DCP/bcca5")

#Tmin
Tmin <- nc_open("Extraction_tasmin.nc") # opens netcdf file example.nc as R object
print (Tmin) # show information about the structure of the data
#Tmax
Tmax <- nc_open("Extraction_tasmax.nc") # opens netcdf file example.nc as R object
print (Tmax) # show information about the structure of the data

lat = ncvar_get(Tmin, "latitude")
lon = ncvar_get(Tmin,"longitude")
time = ncvar_get(Tmin,"time") 
projection= ncvar_get(Tmin,"projection")

tmax = ncvar_get(Tmax,varid="tasmax")
tmin = ncvar_get(Tmin,varid="tasmin")

#remove ncdf
nc_close(Tmin)
nc_close(Tmax)

leap.year=function(year){
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

# Define time: 1950Jan through 2099Dec
years= 1950:2099
J= 1:365
time.mat= expand.grid(J, years)
#add days in leap years
leap.years=years[leap.year(years)==TRUE]
leap.days= expand.grid(366, leap.years)
time.mat= rbind(time.mat, leap.days)
time.mat= time.mat[order(time.mat[,2], time.mat[,1]),]
times= paste(time.mat[,2], time.mat[,1], sep="")

#--------------------------------------------

#RESTRICT ELEVATION  
#DEM
elev= raster::getData('worldclim', var='alt', res=5)
#crop
e <- extent(min(lon)-360, max(lon)-360, min(lat), max(lat))
co.elev= crop(elev, e)

#create grid of elevations  
pts=  expand.grid(lon-360, rev(lat) )
#add indices
pts.ind= expand.grid(1:length(lon), length(lat):1 )
pts= cbind(pts, pts.ind)
names(pts)= c("lon","lat","lon.ind", "lat.ind")
pts$ind= 1:nrow(pts)

pts$elev= raster::extract(co.elev, pts[,1:2])
elev.mat= matrix(pts$elev, nrow=length(lat), ncol=length(lon), byrow=TRUE)

filled.contour(lon,lat, t(elev.mat[length(lat):1,]), color= terrain.colors, asp=1)  #[length(lat):1,]
plot(co.elev)

#Subset to elevation range
elev.mat.sel= elev.mat
elev.mat.sel[elev.mat.sel<1200]= NA
elev.mat.sel[elev.mat.sel>3200]= NA
filled.contour(lon,lat, t(elev.mat.sel[length(lat):1,]), color= terrain.colors, asp=1)  

#=====================================
#calculate initial lambdas (1950) in all cells at reasonable elevations

#ANALYSIS
pts$airpr= sapply(pts$elev, FUN=airpressure_elev) 

#---------------------------------------------
#Load microclimate data

#path for microclimate data
mpath<- '/Volumes/GoogleDrive/My\ Drive/Buckley/Work/Butterflies/Plasticity/ClimateData/microclim/COextract/'

month=6
Tsurf.sun6= read.csv(paste(mpath,"Tsurf_sun_",month,".csv",sep="")) 
Tsurf.shade6= read.csv(paste(mpath,"Tsurf_shade_",month,".csv",sep="")) 
solar6= read.csv(paste(mpath,"solar_",month,".csv",sep="")) 
zenith6= read.csv(paste(mpath,"zenith_",month,".csv",sep="")) 
wind6= read.csv(paste(mpath,"wind_",month,".csv",sep="")) 

month=7
Tsurf.sun7= read.csv(paste(mpath,"Tsurf_sun_",month,".csv",sep="")) 
Tsurf.shade7= read.csv(paste(mpath,"Tsurf_shade_",month,".csv",sep="")) 
solar7= read.csv(paste(mpath,"solar_",month,".csv",sep="")) 
zenith7= read.csv(paste(mpath,"zenith_",month,".csv",sep="")) 
wind7= read.csv(paste(mpath,"wind_",month,".csv",sep="")) 

month=8
Tsurf.sun8= read.csv(paste(mpath,"Tsurf_sun_",month,".csv",sep="")) 
Tsurf.shade8= read.csv(paste(mpath,"Tsurf_shade_",month,".csv",sep="")) 
solar8= read.csv(paste(mpath,"solar_",month,".csv",sep="")) 
zenith8= read.csv(paste(mpath,"zenith_",month,".csv",sep="")) 
wind8= read.csv(paste(mpath,"wind_",month,".csv",sep="")) 

#Load soil temp data for development
Tsurf.shade3= read.csv(paste(mpath,"Tsurf_shade_3.csv",sep=""))
Tsurf.shade4= read.csv(paste(mpath,"Tsurf_shade_4.csv",sep=""))
Tsurf.shade5= read.csv(paste(mpath,"Tsurf_shade_5.csv",sep=""))
Tsurf.shade6= read.csv(paste(mpath,"Tsurf_shade_6.csv",sep=""))
Tsurf.shade9= read.csv(paste(mpath,"Tsurf_shade_9.csv",sep=""))

#-----------------------------------------
#Assemble data

#elevation subset
grid.sel= which(pts$elev>1200 & pts$elev<3200)
pts.sel= pts[grid.sel,1:7]

#Other microclimate 
#Max and min temp in shade
Ts_sh_min= rbind(apply(Tsurf.shade3[pts.sel[,"ind"],3:26], 1, min),apply(Tsurf.shade4[pts.sel[,"ind"],3:26], 1, min),apply(Tsurf.shade5[pts.sel[,"ind"],3:26], 1, min),apply(Tsurf.shade6[pts.sel[,"ind"],3:26], 1, min),apply(Tsurf.shade7[pts.sel[,"ind"],3:26], 1, min),apply(Tsurf.shade8[pts.sel[,"ind"],3:26], 1, min),apply(Tsurf.shade9[pts.sel[,"ind"],3:26], 1, min) ) 

Ts_sh_max= rbind(apply(Tsurf.shade3[pts.sel[,"ind"],3:26], 1, max),apply(Tsurf.shade4[pts.sel[,"ind"],3:26], 1, max),apply(Tsurf.shade5[pts.sel[,"ind"],3:26], 1, max),apply(Tsurf.shade6[pts.sel[,"ind"],3:26], 1, max),apply(Tsurf.shade7[pts.sel[,"ind"],3:26], 1, max),apply(Tsurf.shade8[pts.sel[,"ind"],3:26], 1, max),apply(Tsurf.shade9[pts.sel[,"ind"],3:26], 1, max) ) 


#Other microclimate 
solar=  abind(solar6[pts.sel[,"ind"],3:26],solar7[pts.sel[,"ind"],3:26], solar8[pts.sel[,"ind"],3:26], along=3)

Ts_sun= abind(Tsurf.sun6[pts.sel[,"ind"],3:26],Tsurf.sun7[pts.sel[,"ind"],3:26], Tsurf.sun8[pts.sel[,"ind"],3:26], along=3)
Ts_sh= abind(Tsurf.shade6[pts.sel[,"ind"],3:26],Tsurf.shade7[pts.sel[,"ind"],3:26], Tsurf.shade8[pts.sel[,"ind"],3:26], along=3) 
wind= abind(wind6[pts.sel[,"ind"],3:26], wind7[pts.sel[,"ind"],3:26], wind8[pts.sel[,"ind"],3:26], along=3) 
zenith= abind(zenith6[pts.sel[,"ind"],3:26],zenith7[pts.sel[,"ind"],3:26], zenith8[pts.sel[,"ind"],3:26], along=3) 
zenith[zenith>=80]=80 #set zenith position below horizon to psi=80degrees

### SAVE DATA
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/Data/")
write.csv(pts.sel, "pts.sel.csv")
write.csv(Ts_sh_min, "Ts_sh_min.csv")
write.csv(Ts_sh_max, "Ts_sh_max.csv")

#save objects
saveRDS(solar, "solar.rds")
saveRDS(Ts_sun, "Ts_sun.rds")
saveRDS(Ts_sh, "Ts_sh.rds")
saveRDS(wind, "wind.rds")
saveRDS(zenith, "zenith.rds")



