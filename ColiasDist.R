#========================================
# ISSUES
#Run across years
#Uses constant microclimate

#CHOOSE PROJECTION
proj.k=2 #1: bcc-csm1-1.1.rcp60, 2: ccsm4.1.rcp60, 3: gfdl-cm3.1.rcp60
# 1 has jump and 3 dips in recent
projs=c("bcc-csm","ccsm4","gfdl")
#========================================

memory.size(max=TRUE)
memory.limit(size = 4095)

library(zoo)
#library( fields)
#library( evd)
#library( evdbayes)
#library( ismev) 
library(chron) #convert dates
#library(gdata)
#library(maptools)
#library(epicalc) 
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
#library(doMC)
#registerDoMC(cores=4)

library(reshape2)
library(ggplot2)
library(grid)
library(data.table)

fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\"

#source biophysical model and other functions
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Plasticity\\Analysis\\")
#source("ColiasBiophysMod_20July2015_wShade.R")
#source("ColiasBiophysMod_20July2015.R")
source("ColiasBiophysMod_2Mar2016.R")

#other functions
source("MiscFunctions.R")
source("GDDfunction.R")

#source dtr function
##source("DTRfunction.R")
#source zenith angle calculating function
##source("ZenithAngleFunction.R")
#load radiation functions
source("RadiationModel_10Dec2013.R")

#microclimate model
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Evolution\\MicroclimateModel\\")
source("soil_temp_function_14Aug2014_wShade.R") #CHANGED FROM 27 MAY
#source("soil_temp_function_noextrap_20May2014.R")
source('air_temp_at_height_z_function.R') #calculate air temp at some height z

#define geometric mean
geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

## Register Cluster
cores<-8
cl <- makeCluster(cores)
registerDoParallel(cl)

#--------------------------------------
#Read DCP data
# http://gdo-dcp.ucllnl.org/downscaled_cmip_projections/dcpInterface.html

setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Plasticity\\ClimateData\\DCP\\bcca5")

#Tmin
Tmin <- open.ncdf("Extraction_tasmin.nc") # opens netcdf file example.nc as R object
print (Tmin) # show information about the structure of the data
#Tmax
Tmax <- open.ncdf("Extraction_tasmax.nc") # opens netcdf file example.nc as R object
print (Tmax) # show information about the structure of the data

lat = get.var.ncdf(Tmin, "latitude")
lon = get.var.ncdf(Tmin,"longitude")
time = get.var.ncdf(Tmin,"time") 
projection= get.var.ncdf(Tmin,"projection")

tmax = get.var.ncdf(Tmax,varid="tasmax")
tmin = get.var.ncdf(Tmin,varid="tasmin")

#remove ncdf
rm(Tmin)
rm(Tmax)

##projections
#bcc-csm1-1.1.rcp60
#ccsm4.1.rcp60
#gfdl-cm3.1.rcp60

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

#PLOT
#filled.contour(lon,lat,tmax[,,150,1], color= terrain.colors, asp=1)

#--------------------------------------------

#RESTRICT ELEVATION  
#DEM
library(raster)

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

#-------------------------------------------------
#PARAMETERS

#Demographic parameters
#Kingsolver 1983 Ecology 64
OviRate=0.73; # Ovipositing rates: 0.73 eggs/min (Stanton 1980) 
MaxEggs=700; # Max egg production: 700 eggs (Tabashnik 1980)
PropFlight= 0.5; # Females spend 50# of available activity time for oviposition-related
# Watt 1979 Oecologia
SurvDaily=0.6; # Daily loss rate for Colias p. eriphyle at Crested Butte, female values
#Hayes 1981, Data for Colias alexandra at RMBL
SurvMat=0.014; #1.4# survival to maturity

#read/make species data
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Evolution\\Data\\")
#SpecDat<-read.csv("SpecDat_23July2013.csv");
solar.abs= 0.65 #seq(0.4,0.7,0.01) #seq(0.4,0.7,0.001)
SpecDat= as.data.frame(solar.abs)
#SpecDat$Species="Cmeadii"
SpecDat$d=0.36
SpecDat$fur.thickness=1.46
SpecDat= as.matrix(SpecDat)
# Solar absorptivity, proportion
# D- Thoractic diameter, cm
# Fur thickness, mm

#Fur thickness from Kingsolver (1983)
#Eriphyle Montrose 0.82mm, d=3.3mm, 1.7km
#Eriphyle Skyland 1.08, d=3.5mm, 2.7km
#Meadii Mesa sco 1.46mm, d=3.6mm, 3.6km
FT= 1.46 #c(0.01, 0.82, 1.46, 2)
#FUR THICKNESS: 0, eriphyle olathe, meadii mesa seco, 2

days<-c(31,28,31,30,31,30,31,31,30,31,30,31) #days in months

#absorptivity
abs1=seq(0.4,0.7,0.05)

#------------------------
#PERFORMANCE FUNCTIONS
#function for flight probability
fl.ph<- function(x) 1 * exp(-0.5*(abs(x-33.5)/5)^3.5)

#function for egg viability
#egg.viab<-function(x) ifelse(x<40, egg.viab<-1, egg.viab<- exp(-(x-40)))
#Set P_survival=0.868 at 45 (Kingsolver dissertation), log(0.868)=-5/t, t=35.32

egg.viab<-function(x) ifelse(x<40, egg.viab<-1, egg.viab<- exp(-(x-40)/35.32))
#plot(30:50, egg.viab(30:50))

#estimate flight time for Meadii using data from Ward 1977
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\Butterflies\\Evolution\\Data\\")
bpop= read.csv("ButterflyPop_Ward1977.csv")
bpop= subset(bpop, bpop$Loc=="CumbPass")
bpop.tot= sum(bpop$Total)
bpop$prop= bpop$Total/bpop.tot

#Calculate weighted mean and sd for flight date
flightday.mean=wt.mean(bpop$Date, bpop$Total)
flightday.sd=wt.sd(bpop$Date, bpop$Total)

air_temp_at_height_z.mat<- function(T_mat, z_0, z_r, z){
  T_r= T_mat[1]
  T_s= T_mat[2]
  T_z<-(T_r-T_s)*log((z+z_0)/z_0+1)/log((z_r+z_0)/z_0+1)+T_s ##this is exactly eqn (19) of the notes
  return(T_z)
}

#---------------------------------------------
#Load microclimate data

#path for microclimate data
mpath<- 'C:/Users/Buckley/Google Drive/Buckley/Work/Butterflies/Plasticity/ClimateData/microclim/COextract/'

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

#---------------------------------------
#Temps at plant height
z_0_1=0.02

mo= c(rep(1,31),rep(2,28),rep(3,31), rep(4,30),rep(5,31),rep(6,30),rep(7,31),rep(8,31),rep(9,30) )
# 60:273, march 1 to sep 30
#=========================================================================

#SET PARAMETERS

### Across sites
# SET UP CLIMATE DATA
spec.k=1
ft.k=1
D=SpecDat[1,"d"]; delta= SpecDat[1,"fur.thickness"]

#---------------------------
#Make array to store data
#Lambda= survival*fecundity
#Matrix of lambdas
#dims: stats, year, Lambda 
Lambda<-array(NA, dim=c(length(years),nrow(pts.sel),length(seq(0.4,0.7,0.05)),3,4)) #Last dimension is Lambda, FAT,Egg Viability
dimnames(Lambda)[[1]]<-years
#dimnames(Lambda_20cm)[[2]]<-stats
#dimnames(Lambda)[[3]]="abs"
dimnames(Lambda)[[4]]<-c("gen1","gen2","gen3")

#Matrix for pupual temps
pup.temps<-array(NA, dim=c(12, length(years),nrow(pts.sel),3)) #3 generations 
#Add names
dimnames(pup.temps)[[1]]= c("stat","yr","gen","Jlarv", "Jpup","Jadult","Tlarv","Tpup","Tad","Tlarv_fixed","Tpup_fixed","Tad_fixed") 

#Te
dayk= 152:243
Te.mat.all= array(data=NA, dim=c(nrow(pts.sel), length(6:20),7, length(dayk) ) )

aseq= seq(0.4,0.7,0.05)

#===================================================================
#LOOP YEARS

for(yr.k in 1:length(years) ){
print(yr.k)
  
inds= which(time.mat[,2]==years[yr.k])

tmax.yr= tmax[, , inds, proj.k]  
tmin.yr= tmin[, , inds, proj.k] 

#============================================================================
#CALCULATE DEVELOPMENT TIMING AND TEMPS

#Temps at plant height
lon.inds= pts.sel[, "lon.ind"]
lat.inds= pts.sel[, "lat.ind"]
lonlat.inds <- cbind(lon.inds,lat.inds)

Ta_plant_min<- foreach(d=60:273, .combine='cbind')  %do% {
  T_mat= tmin.yr[, ,d]
  T_mat= cbind(T_mat[lonlat.inds], Ts_sh_min[mo[d]-2, ])
  apply(T_mat,MARGIN=1, FUN=air_temp_at_height_z.mat, z_0=z_0_1, z_r=2, z=z_0_1)}

Ta_plant_max<- foreach(d=60:273, .combine='cbind')  %do% {
  T_mat= tmax.yr[, ,d]
  T_mat= cbind(T_mat[lonlat.inds], Ts_sh_max[mo[d]-2, ])
  apply(T_mat,MARGIN=1, FUN=air_temp_at_height_z.mat, z_0=z_0_1, z_r=2, z=z_0_1)}

#transpose
Ta_plant_min= t(Ta_plant_min)
Ta_plant_max= t(Ta_plant_max)

##subset post snowment ## ELEVATION REGRESSION ACROSS WESTERN CLIMATE DATA CENTER?
#dat2= subset(dat2, dat2$Julian>= Jmelt)

#----------------------------
# ESTIMATE DEVELOPMENTAL TIMING

    DevZeros= c(9.2176, 11.5, 9.7) #4th and 5th, larval, pupal
    GddReqs= c(117.06, 270.39 ,101.9) 
    
    #Calc GDDs
    GDDs_45<- foreach(d=(60:273)-59, .combine='cbind')  %do% {
      T_mat= cbind(Ta_plant_min[d,], Ta_plant_max[d,])
      apply(T_mat,MARGIN=1, FUN=degree.days.mat, LDT=DevZeros[1])}
      
    GDDs_l<- foreach(d=(60:273)-59, .combine='cbind')  %do% {
      T_mat= cbind(Ta_plant_min[d,], Ta_plant_max[d,])
      apply(T_mat,MARGIN=1, FUN=degree.days.mat, LDT=DevZeros[2])}
    
    GDDs_p<- foreach(d=(60:273)-59, .combine='cbind')  %do% {
      T_mat= cbind(Ta_plant_min[d,], Ta_plant_max[d,])
      apply(T_mat,MARGIN=1, FUN=degree.days.mat, LDT=DevZeros[3])}
    
    GDDs_45= t(GDDs_45)
    GDDs_l= t(GDDs_l)
    GDDs_p= t(GDDs_p)

    #-----------------------------------
    #Calculate development timing
    Js=60:273
    
for(cell.k in 1:nrow(pts.sel)){  #  for(cell.k in 1:nrow(pts.sel))
    
Ta_plant_mean= rowMeans(cbind(Ta_plant_min[, cell.k], Ta_plant_max[, cell.k]))

for(gen.k in 1:3){
      
      #Assume 7 days from eclosion to eggs laid
      #Hatching ~5days (~70hrs, based on Heidi's heat shock data, Jessica's development data doesn't seem to have hatching time)
      Jlarv= ifelse(gen.k>1, Jadult+12, Js[which.max(GDDs_45[,cell.k]>0)] )  
      if(Jlarv>max(Js)) Jlarv=max(Js)
      
      ##TO PUPATION
      check=ifelse(gen.k==1, gdds<-GDDs_45[,cell.k],gdds<-GDDs_l[,cell.k])
      
      Jpup= Jlarv + which.max(cumsum(gdds[which(Js==Jlarv):length(Js)])> ifelse(gen.k==1, GddReqs[1],GddReqs[2])  )
      if(Jpup>max(Js) | length(Jpup)==0) Jpup=max(Js) 
      
      #PUPATION
      gdds<-GDDs_p[,cell.k]
      Jadult= Jpup + which.max(cumsum(gdds[which(Js==Jpup):length(Js)])> GddReqs[3])
      if(Jadult>max(Js) | length(Jadult)==0) Jadult=max(Js)
      
      #----------------------
      #Calculate temps
      Tlarv= mean(Ta_plant_mean[Js %in% Jlarv:Jpup], na.rm=TRUE)
      Tpup= mean(Ta_plant_mean[Js %in% Jpup:Jadult], na.rm=TRUE)
      Tad= mean(Ta_plant_mean[Js %in% Jadult:(Jadult+7)], na.rm=TRUE)
      ### ADULT TEMP IS AIR TEMP
      #Check if more than 5 NAs
          
      #Write data in array
      pup.temps[3:9,yr.k, cell.k, gen.k]=c(gen.k,Jlarv,Jpup,Jadult,Tlarv,Tpup,Tad)
      
    } #end loop generation
    
} #end loop cells

#=========================================================================
#Run Te calculations

#Set up climate data
hr.dec= 1:24

tmax.yr= tmax[, , inds[152:243], proj.k]  
tmin.yr= tmin[, , inds[152:243], proj.k]  

clim= as.data.frame(152:243)
names(clim)="J"
clim$month=NA
clim[which(clim$J<182),"month"]=6
clim[which(clim$J>181 & clim$J<213),"month"]=7
clim[which(clim$J>212),"month"]=8

#estimate daylength
Trise.set= suncalc(clim$J, Lat = pts.sel[1,"lat"], Long = pts.sel[1,"lon"], UTC = FALSE)
clim$set= Trise.set$sunset
clim$rise= Trise.set$sunrise

#TEMP
Thr<- foreach(cell.k=1:length(grid.sel) ) %do% {
  t(apply( cbind( tmax.yr[pts.sel[cell.k, "lon.ind"], pts.sel[cell.k, "lat.ind"],], tmin.yr[pts.sel[cell.k, "lon.ind"], pts.sel[cell.k, "lat.ind"],], clim[,c("rise","set")]) , FUN=Thours.mat, MARGIN=1))}
#Collapse list into array
Thr= array(unlist(Thr), dim = c(nrow(Thr[[1]]), ncol(Thr[[1]]), length(Thr)))

#RADIATION
Rhr= foreach(cell.k=1:length(grid.sel)) %do% {
  t(apply(rbind(matrix(solar6[pts.sel[cell.k,"ind"],3:26],nrow = 30,ncol = 24, byrow=TRUE), matrix(solar7[pts.sel[cell.k,"ind"],3:26],nrow = 31,ncol = 24, byrow=TRUE), matrix(solar8[pts.sel[cell.k,"ind"],3:26],nrow = 31,ncol = 24, byrow=TRUE) ), FUN=Rad.mat, MARGIN=1)) } 
#Collapse list into array
Rhr= array(unlist(Rhr), dim = c(nrow(Rhr[[1]]), ncol(Rhr[[1]]), length(Rhr)))
#columns 1:24 are direct, 25:28 are diffuse
#reflected is direct * albedo of 0.7

#Other microclimate 
Ts_sun= abind(Tsurf.sun6[pts.sel[,"ind"],3:26],Tsurf.sun7[pts.sel[,"ind"],3:26], Tsurf.sun8[pts.sel[,"ind"],3:26], along=3)
Ts_sh= abind(Tsurf.shade6[pts.sel[,"ind"],3:26],Tsurf.shade7[pts.sel[,"ind"],3:26], Tsurf.shade8[pts.sel[,"ind"],3:26], along=3) 
wind= abind(wind6[pts.sel[,"ind"],3:26], wind7[pts.sel[,"ind"],3:26], wind8[pts.sel[,"ind"],3:26], along=3) 
zenith= abind(zenith6[pts.sel[,"ind"],3:26],zenith7[pts.sel[,"ind"],3:26], zenith8[pts.sel[,"ind"],3:26], along=3) 
zenith[zenith>=80]=80 #set zenith position below horizon to psi=80degrees

#--------------------------------------------------
#Calculate Te #cell.k: length(grid.sel)

for(hr.k in 1:15 ){
  
  Thr.d= Thr[,hr.k,]
  Ts_sun.d= Ts_sun[,hr.k,]
  Ts_sh.d= Ts_sh[,hr.k,]
  wind.d= wind[,hr.k,]
  Rhr.d= Rhr[,c(hr.k,hr.k+24),]
  zenith.d= zenith[,hr.k,]
  
  #combine data
  Te.dat=abind(t(Thr.d), Ts_sun.d[,(clim$month-5)], Ts_sh.d[,(clim$month-5)], wind.d[,(clim$month-5)], t(Rhr.d[,1,]), t(Rhr.d[,2,]),zenith.d[,(clim$month-5)],along=3)
  
  for(a.k in 1:7){
  
  a= aseq[a.k]
  Te.mat.all[,hr.k,a.k,]<-  apply(Te.dat,MARGIN=c(1,2), FUN=biophys.var_sh.mat, D, delta, a)  
  
  } # end loop across absorptivities
} # end loop across hours
 
#=======================================
#DEMOGRAPHY

#Flight probability
FlightProb<- foreach(a= 1:dim(Te.mat.all)[[3]] )  %:% foreach(hr=1:15) %:% foreach(d=1:dim(Te.mat.all)[[4]], .combine='cbind') %do% {  out=sapply(Te.mat.all[,hr, a,d], FUN=fl.ph) }
#Back to array
FlightProb.all= array(unlist(FlightProb), dim = c(nrow(FlightProb[[1]][[1]]), ncol(FlightProb[[1]][[1]]), length(FlightProb[[1]]), length(FlightProb)))

#Egg viability
EggViab<- foreach(a= 1:dim(Te.mat.all)[[3]] )  %:% foreach(hr=1:15) %:% foreach(d=1:dim(Te.mat.all)[[4]], .combine='cbind') %do% {  out=sapply(Te.mat.all[,hr, a,d], FUN=egg.viab) }
#Back to array        
EggViab.all= array(unlist(EggViab), dim = c(nrow(EggViab[[1]][[1]]), ncol(EggViab[[1]][[1]]), length(EggViab[[1]]), length(EggViab)))

#-------------------------------------------
#DEMOGRAPHY

EV1=matrix(NA,2, nrow(SpecDat))

  for(gen.k in 1:3 ){ #loop generation
  
    #find cells that can complete generation
    Jfls=pup.temps["Jadult",yr.k,, gen.k] 
    cell.inds= which(Jfls<244 & Jfls>151)
    
    for(cell.k in cell.inds ){ #loop cells
      # if(cell.k/100==round(cell.k/100)) print(cell.k)
      
      for(abs.k in 1:dim(Te.mat.all)[3] ){ #loop absorptivity
        
      Te.mat= Te.mat.all[cell.k,,abs.k,]
      
      #get flight dates
      Jfl=pup.temps["Jadult",yr.k, cell.k, gen.k]   
   
    #average over hours
    FAT= rowSums(FlightProb.all[cell.k,,,abs.k], na.rm=TRUE)
    EggViab= apply(EggViab.all[cell.k,,,abs.k], FUN=geo_mean, MARGIN=1) #Egg viability GEOMETRIC MEAN ACROSS HOURS
    Temps= colMeans(Te.mat[,], na.rm=TRUE)
    
    ##CALCULATE EGG VIABILITY OVER 5 DAY PERIOD (GEOMETRIC MEAN ACROSS HOURS)
    #sample flight day from truncated normal distribution
    Nind=1000 #changed from 100
    f.low= max(Jfl-7,min(clim$J)+2)
    f.up= min(Jfl+7,max(clim$J)-2)
    
    flightday= round(rtruncnorm(Nind, a=f.low, b=f.up, mean = Jfl, sd = 2) )
    f.ind= match(flightday, clim$J)
    #if NA day, use mean
    f.ind[is.na(f.ind)]<-match(Jfl, clim$J)
    
    #calculate geometric mean of egg viability within flight period
    ev.ind=sapply(f.ind, function(x)  geo_mean(EggViab[(x-2):(x+2)]) )
    #AVERAGE FAT OVER DAYS
    FAT.ind= sapply(f.ind, function(x)  mean(FAT[(x-2):(x+2)], na.rm=TRUE) )
    #AVERAGE TEMP
    T.ind= sapply(f.ind, function(x)  mean(Temps[(x-2):(x+2)], na.rm=TRUE) )
    
    Eggs.ind= 60*PropFlight*OviRate*FAT.ind * ev.ind #account for Egg viability
    Eggs.ind_noViab= 60*PropFlight*OviRate*FAT.ind
    
    #Means across individuals
    Eggs= mean(Eggs.ind)
    Eggs_noViab= mean(Eggs.ind_noViab)
    EV1[1,spec.k]= mean(FAT.ind)
    EV1[2,spec.k]= mean(ev.ind)
    
    if(!is.nan(Eggs)){
      MaxDay=5
      Lambda1=0
      for(day in 1:MaxDay){
        Eggs1= min(Eggs, MaxEggs-Eggs*(day-1))  ###LIMIT MAX NUMBER EGGS
        if(Eggs1<0) Eggs1=0
        Lambda1= Lambda1+ SurvMat * SurvDaily^day *Eggs1;                        
      }#end loop days
      
      Lambda[yr.k, cell.k, abs.k, gen.k, ]= c(Lambda1, mean(FAT.ind), mean(ev.ind), mean(T.ind, na.rm=T) ) 
     
    } #Check Eggs

      } #end loop absortivity
    } #end loop cells    
  } #end loop generation
  
#Save every 30 years
if(yr.k/30== round(yr.k/30)){
#SAVE OBJECT
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")

filename= paste("lambda1_",projs[proj.k],".rds",sep="")
saveRDS(Lambda, filename)
#Lambda1 <- readRDS("mymodel.rds")

#Write out pupal temps
saveRDS(pup.temps, paste("PupTemps_",projs[proj.k],".rds",sep="") )
} #end save

} #end loop years

#SAVE OBJECT
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")

filename= paste("lambda1_",projs[proj.k],".rds",sep="")
saveRDS(Lambda, filename)
#Lambda1 <- readRDS("mymodel.rds")

#Write out pupal temps
saveRDS(pup.temps, paste("PupTemps_",projs[proj.k],".rds",sep="") )

#write out points
write.csv(pts.sel, paste("COpoints_",projs[proj.k],".rds",sep="") )

#========================================
#PLOT FITNESS CURVES
plot(seq(0.4,0.7,0.05), Lambda[1, 10, , 1, 1])
plot(seq(0.4,0.7,0.05), Lambda[1, 10, , 1, 2])

#==========================================
#PLOT TEMPERATURE TRENDS

#filled.contour(lon,lat, t(elev.mat.sel[length(lat):1,]), color= terrain.colors, asp=1) 

#plot time series
plot(tmax[15,15,,proj.k], type="l")
plot(tmin[15,15,,proj.k], type="l")

projs=c("bcc-csm","ccsm4","gfdl")

#loop projections
for(proj.k in 1:3){

#annual averaging
tmin.ann= foreach(cell.k=1:nrow(pts.sel), .combine='cbind') %do% aggregate(tmin[pts.sel[cell.k, "lon.ind"], pts.sel[cell.k, "lat.ind"],,proj.k], list(time.mat[,2]), FUN=base::mean)[,2]

tmax.ann= foreach(cell.k=1:nrow(pts.sel), .combine='cbind') %do% aggregate(tmax[pts.sel[cell.k, "lon.ind"], pts.sel[cell.k, "lat.ind"],,proj.k], list(time.mat[,2]), FUN=base::mean)[,2]

setwd(paste(fdir,"figures\\",sep="") )
pdf(paste("Tannual_acrossSites_",projs[proj.k],".pdf", sep=""), height = 6, width = 6)
plot(years, tmin.ann[,4], type="l", ylab="T annual mean", ylim=range(-7,20))
points(years, tmax.ann[,4], type="l", lty="dashed")
legend("topleft", c("Tmax","Tmin"),lty=c("dashed", "solid"))
dev.off()

#add names
colnames(tmin.ann)= pts.sel$elev
colnames(tmax.ann)= pts.sel$elev

#melt data for ggplot
tmin2 <- melt(tmin.ann, variable.name ='elev', value.name='temp')
tmin2[,1]= years[tmin2[,1]]
colnames(tmin2)[1:2]=c("yr","elev")
#cut by elevation
tmin2$ecut= cut(tmin2$elev, breaks=3)

t1= ggplot(tmin2, aes(x=yr, y=temp, group=elev, color=elev) ) +geom_smooth(method=loess,se=FALSE)  +facet_grid(~ecut)

tmax2 <- melt(tmax.ann, variable.name ='elev', value.name='temp')
tmax2[,1]= years[tmax2[,1]]
colnames(tmax2)[1:2]=c("yr","elev")
#cut by elevation
tmax2$ecut= cut(tmax2$elev, breaks=3)

t2= ggplot(tmax2, aes(x=yr, y=temp, group=elev, color=elev) ) +geom_smooth(method=loess,se=FALSE) +facet_grid(~ecut)

#plot temp across all sites
setwd(paste(fdir,"figures\\",sep="") )
pdf(paste("Tannual_allSites_",projs[proj.k],".pdf", sep=""), height = 12, width = 12)

grid.newpage()
pushViewport(viewport(layout=grid.layout(2,1)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(t1,vp=vplayout(1,1))
print(t2,vp=vplayout(2,1))

dev.off()

} #end projection loop

#--------------     
#summer averaging
summ.inds= which(time.mat[,1]>165 & time.mat[,1]<259)

tmin.ann= foreach(cell.k=1:nrow(pts.sel), .combine='cbind') %do% aggregate(tmin[pts.sel[cell.k, "lon.ind"], pts.sel[cell.k, "lat.ind"],summ.inds, proj.k], list(time.mat[summ.inds,2]), FUN=base::mean)[,2]
tmax.ann= foreach(cell.k=1:nrow(pts.sel), .combine='cbind') %do% aggregate(tmax[pts.sel[cell.k, "lon.ind"], pts.sel[cell.k, "lat.ind"],summ.inds, proj.k], list(time.mat[summ.inds,2]), FUN=base::mean)[,2]

#add names
colnames(tmin.ann)= pts.sel$elev
colnames(tmax.ann)= pts.sel$elev

#melt data for ggplot
tmin2 <- melt(tmin.ann, variable.name ='elev', value.name='temp')
tmin2[,1]= years[tmin2[,1]]
colnames(tmin2)[1:2]=c("yr","elev")
ggplot(tmin2, aes(x=yr, y=temp, group=elev, color=elev) ) +geom_smooth(method=loess,se=FALSE)
#ggplot(tmin2, aes(x=yr, y=temp, group=elev, color=elev) ) +geom_line()

tmax2 <- melt(tmax.ann, variable.name ='elev', value.name='temp')
tmax2[,1]= years[tmax2[,1]]
colnames(tmax2)[1:2]=c("yr","elev")
ggplot(tmax2, aes(x=yr, y=temp, group=elev, color=elev) ) +geom_smooth(method=loess,se=FALSE)
#ggplot(tmax2, aes(x=yr, y=temp, group=elev, color=elev) ) +geom_line()

#============================================================
#Load and check pupal tempertures and lambdas

#Read lambdas and pupal temps
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\3gen_rds")

Lambda_old <- readRDS( paste("lambda1_",projs[proj.k],".rds", sep="") )
#last dim: lambda, FAT, egg viability, temp

pup.temps <- readRDS( paste("PupTemps_",projs[proj.k],".rds", sep="") )
#first dim:  stat          yr         gen       Jlarv        Jpup      Jadult       Tlarv  Tpup         Tad Tlarv_fixed  Tpup_fixed   Tad_fixed 

summary(pup.temps[,1,,3])

na_count <-function (x) sapply(x, function(y) sum(is.na(y)))

#No NAs in pupal temps
nas= na_count(pup.temps[9,,,3])
sum(nas)

# NAs in lambdas
nas= na_count(Lambda[9,,3,1,4])
sum(nas)

#------------
#pts.sel # 841   8
#Lambda[yr.k, , , gen.k, ] #150 841   7   3   4
#pup.temps[yr.k, , , gen.k, ] # 12 150 841   3

lamb= cbind(pts.sel,t(Lambda[1:20, , 3, 2, 1]) )

pup.temps[6, 1:20, , 2]

Lambda_old[1, , 3, 2, 1]
Lambda[1, , 3, 2, 1]

#--------------------
#PLOTS

lper1= cbind(pts.sel, Ta_plant_min[1,])
lper1= cbind(pts.sel, GDDs_l[1,])

lper1= cbind(pts.sel, pup.temps[9,yr.k,,2]) 
names(lper1)[8]="var"
quilt.plot(lper1$lon,lper1$lat, lper1$var)

plot(lper1$lon, lper1$var)



