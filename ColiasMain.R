
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
library(truncnorm)
#========================================

#CHOOSE PROJECTION
proj.k=2 #1: bcc-csm1-1.1.rcp60, 2: ccsm4.1.rcp60, 3: gfdl-cm3.1.rcp60
projs=c("bcc-csm","ccsm4","gfdl")

#LOAD PARAMETERS
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
solar.abs= 0.65 # Solar absorptivity, proportion
SpecDat= as.data.frame(solar.abs)
SpecDat$d=0.36 # d- Thoractic diameter, cm
SpecDat$fur.thickness=1.46 # Fur thickness, mm
SpecDat= as.matrix(SpecDat)

#Fur thickness from Kingsolver (1983)
#Eriphyle Montrose 0.82mm, d=3.3mm, 1.7km
#Eriphyle Skyland 1.08, d=3.5mm, 2.7km
#Meadii Mesa sco 1.46mm, d=3.6mm, 3.6km
FT= 1.46 #c(0.01, 0.82, 1.46, 2)
#FUR THICKNESS: 0, eriphyle olathe, meadii mesa seco, 2

days<-c(31,28,31,30,31,30,31,31,30,31,30,31) #days in months

#absorptivity
abs1=seq(0.4,0.7,0.05)

#flight time for Meadii using data from Ward 1977
#weighted mean and sd for flight date
flightday.mean=210.95
flightday.sd=7.27

#Temps at plant height
z_0_1=0.02

mo= c(rep(1,31),rep(2,28),rep(3,31), rep(4,30),rep(5,31),rep(6,30),rep(7,31),rep(8,31),rep(9,30) )
# 60:273, march 1 to sep 30

# SET UP CLIMATE DATA
spec.k=1
ft.k=1
D=SpecDat[1,"d"]; delta= SpecDat[1,"fur.thickness"]

#===========================
#LOAD OR DEFINE FUNCTIONS

#specify directory
fdir= ("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/") 
 
#microclimate model
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/Butterflies/Evolution/MicroclimateModel/")
source("soil_temp_function_14Aug2014_wShade.R")

#radiation
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/Butterflies/Plasticity/Analysis/")
source("RadiationModel_10Dec2013.R")

#define geometric mean
geo_mean <- function(data) {
  log_data <- log(data)
  gm <- exp(mean(log_data[is.finite(log_data)]))
  return(gm)
}

#PERFORMANCE FUNCTIONS
#function for flight probability
fl.ph<- function(x) 1 * exp(-0.5*(abs(x-33.5)/5)^3.5)

#function for egg viability
egg.viab<-function(x) ifelse(x<40, egg.viab<-1, egg.viab<- exp(-(x-40)/35.32))


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

# Define time: 1950Jan through 2099Dec
leap.year=function(year){
  return(((year %% 4 == 0) & (year %% 100 != 0)) | (year %% 400 == 0))
}

years= 1950:2099
J= 1:365
time.mat= expand.grid(J, years)
#add days in leap years
leap.years=years[leap.year(years)==TRUE]
leap.days= expand.grid(366, leap.years)
time.mat= rbind(time.mat, leap.days)
time.mat= time.mat[order(time.mat[,2], time.mat[,1]),]
times= paste(time.mat[,2], time.mat[,1], sep="")

#=========================================================================

### LOAD DATA
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/Data/")
pts.sel<- read.csv("pts.sel.csv")
Ts_sh_min<- read.csv("Ts_sh_min.csv")
Ts_sh_max<- read.csv("Ts_sh_max.csv")

#load objects
solar <- readRDS("solar.rds")
Ts_sun <- readRDS("Ts_sun.rds")
Ts_sh <- readRDS("Ts_sh.rds")
wind <- readRDS("wind.rds")
zenith <- readRDS("zenith.rds")

### SET UP DATA STRUCTURES
#Make array to store data
#Lambda= survival*fecundity
#Matrix of lambdas
#dims: stats, year, Lambda 
Lambda<-array(NA, dim=c(length(years),nrow(pts.sel),length(seq(0.4,0.7,0.05)),3,4)) #Last dimension is Lambda, FAT,Egg Viability
dimnames(Lambda)[[1]]<-years
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
    apply(T_mat,MARGIN=1, FUN=air_temp_at_height_z_mat, z_0=z_0_1, z_r=2, z=z_0_1)}
  
  Ta_plant_max<- foreach(d=60:273, .combine='cbind')  %do% {
    T_mat= tmax.yr[, ,d]
    T_mat= cbind(T_mat[lonlat.inds], Ts_sh_max[mo[d]-2, ])
    apply(T_mat,MARGIN=1, FUN=air_temp_at_height_z_mat, z_0=z_0_1, z_r=2, z=z_0_1)}
  
  #transpose
  Ta_plant_min= t(Ta_plant_min)
  Ta_plant_max= t(Ta_plant_max)
  
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
  Trise.set= suncalc(clim$J, Lat = pts.sel[1,"lat"], Long = pts.sel[1,"lon"])
  clim$set= Trise.set$sunset
  clim$rise= Trise.set$sunrise
  
  #TEMP
  Thr<- foreach(cell.k=1:nrow(pts.sel) ) %do% {
    t(apply( cbind( tmax.yr[pts.sel[cell.k, "lon.ind"], pts.sel[cell.k, "lat.ind"],], tmin.yr[pts.sel[cell.k, "lon.ind"], pts.sel[cell.k, "lat.ind"],], clim[,c("rise","set")]) , FUN=Thours.mat, MARGIN=1))}
  #Collapse list into array
  Thr= array(unlist(Thr), dim = c(nrow(Thr[[1]]), ncol(Thr[[1]]), length(Thr)))
  
  #RADIATION
  Rhr= foreach(cell.k=1:nrow(pts.sel)) %do% {
    t(apply(rbind(matrix(solar[cell.k,,1],nrow = 30,ncol = 24, byrow=TRUE), matrix(solar[cell.k,,2],nrow = 31,ncol = 24, byrow=TRUE), matrix(solar[cell.k,,3],nrow = 31,ncol = 24, byrow=TRUE) ), FUN=Rad.mat, MARGIN=1)) } 
  #Collapse list into array
  Rhr= array(unlist(Rhr), dim = c(nrow(Rhr[[1]]), ncol(Rhr[[1]]), length(Rhr)))
  #columns 1:24 are direct, 25:28 are diffuse
  #reflected is direct * albedo of 0.7
  
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
  
} #end loop years

#SAVE OBJECT
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/OUT/")

filename= paste("lambda1_",projs[proj.k],".rds",sep="")
saveRDS(Lambda, filename)
#Lambda1 <- readRDS("mymodel.rds")

#Write out pupal temps
saveRDS(pup.temps, paste("PupTemps_",projs[proj.k],".rds",sep="") )

#write out points
write.csv(pts.sel, paste("COpoints_",projs[proj.k],".rds",sep="") )

#========================================

#EVOLUTIONARY MODEL
N.ind=1000
a= seq(0.4,0.7,0.05)
#for finding a with max fitness
a.fit= as.data.frame(seq(0.4,0.7,0.01))
names(a.fit)="a"

#Read points
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/OUT/") #mac version

pts.sel= read.csv( paste("COpoints.csv", sep="") ) #_",projs[proj.k],"

#Read lambdas and pupal temps
#Lambda[years, sites, abs, gen, metrics: Lambda, FAT,Egg Viability]
Lambda <- readRDS( paste("lambda1_",projs[proj.k],".rds", sep="") )
pup.temps <- readRDS( paste("PupTemps_",projs[proj.k],".rds", sep="") )

#Find years with calculations
counts= rowSums(is.na(pup.temps[6,,,1]))

inds=1:150
years= years[inds]

#==============================================
#Calculate optimal absorptivity

abs.opt= array(NA, dim=c(length(years),nrow(pts.sel), 3))  

for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {
    
    Lambda.yr.gen= Lambda[yr.k, , , gen.k, ]
    
    #Extract temperatures
    Tp= pup.temps["Tpup",yr.k, , gen.k]
    
    #--------------------------
    #Fitness models
    
    if(!all(is.na(Lambda.yr.gen[,,1]))){ #check has data
      
      #Estimate fitness functions across cells
      fit= array(unlist(apply(Lambda.yr.gen[,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2))$coefficients)), dim=c(3, nrow(pts.sel)) )
      #Save model
      fit.mod= apply(Lambda.yr.gen[,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2)) )
      
      #find maxima lambda
      abs.opt[yr.k,,gen.k]= as.vector(array(unlist(sapply(fit.mod, function(x) if(!is.null(x))a.fit$a[which.max(predict.lm(x, a.fit))] )), dim=c(1, nrow(pts.sel)) ) )
      
    } #end check data
  } #end gen loop
} #end year loop

#save optimal Alphas
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/OUT/")

#saveRDS(abs.opt, paste("abs.opt_",projs[proj.k],".rds", sep=""))
abs.opt <- readRDS( paste("abs.opt_",projs[proj.k],".rds", sep="") )

#***************************************
#compute initial AbsMean 
int_elev = 0.4226; slope_elev = 0.06517
Tmid = 20; slope_plast = -0.0083  #if Tmid=22.5, -0.006667;

elev_km= pts.sel$elev/1000
abs.init <- int_elev+ slope_elev*elev_km

## NEED TO CALC ABS.OPT
#initialize with optimum value yrs 1950-1960, across generations
abs.init2 <- rowMeans(colMeans(abs.opt[1:10,, ], na.rm=TRUE))

plot(elev_km, abs.init, ylim=range(0.5, 0.7), type="l")
points(elev_km, abs.init2)

#Use optimal
abs.init<- abs.init2

#-----------------------
#Save values
abs.mean= array(NA, dim=c(length(years),nrow(pts.sel), 3, 5,5))  #dims: yr.k, cell.k, gen.k, scen.k:no plast, plast, only plast, metrics: abssample, absmid, rn, Babsmid, Brn)
abs.mean[1,,1,,2]= abs.init
abs.mean[1,,1,,3]= slope_plast
dimnames(abs.mean)[[5]]= c("abssample", "absmid", "rn", "Babsmid", "Brn") 

lambda.mean= array(NA, dim=c(length(years),nrow(pts.sel), 3, 5)) #dims: yr.k, cell.k, gen.k, scen.k:no plast, plast, only plast)

BetaRN= rep(NA, nrow(pts.sel))
#-------------------------------
scen.mat= rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(1,1,0),c(1,1,1) )
colnames(scen.mat)= c("plast","evol","evolRN"  )

for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {
    
    BetaAbsmid=NA
    
    Lambda.yr.gen= Lambda[yr.k, , , gen.k, ]
    
    #determine those completing generations
    comp.gen= which(pup.temps["Jadult",yr.k,,gen.k]<243)
    nocomp.gen= which(pup.temps["Jadult",yr.k,,gen.k]==243)
    #set those not completing generations to NA
    if(length(nocomp.gen)>0) Lambda.yr.gen[nocomp.gen,,]=NA
    
    #account for NA lambdas
    l.no.na= which(!is.na(Lambda.yr.gen[,1,1]))
    
    if(length(l.no.na)>0){ #CHECK LAMBDA DATA EXISTS
      
      #Extract temperatures
      Tp= pup.temps["Tpup",yr.k,l.no.na, gen.k]
      
      #--------------------------
      #Fitness models
      #Estimate fitness functions across cells
      fit= array(unlist(apply(Lambda.yr.gen[l.no.na,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2))$coefficients)), dim=c(3, nrow(pts.sel)) )
      #Save model
      fit.mod= apply(Lambda.yr.gen[l.no.na,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2)) )
      
      #find maxima lambda
      abs.max= as.vector(array(unlist(sapply(fit.mod, function(x) if(!is.null(x))a.fit$a[which.max(predict.lm(x, a.fit))] )), dim=c(1, length(l.no.na)) ) )
      
      #-------------------------
      # LOOP PLASTICITY SCENARIOS
      for(scen.k in 1:5){ #plast0evol0, plast1evol0, plast0evol1, plast1evol1, plast1evol1rnevol1
        
        if(scen.mat[scen.k,1]==1) rn.mean1= rep(slope_plast, length(l.no.na) )
        if(scen.mat[scen.k,1]==0) rn.mean1= rep(0, length(l.no.na) )
        if(scen.k==5 & gen.k==1) rn.mean1= abs.mean[yr.k,l.no.na,gen.k,scen.k,"rn"]
        if(scen.k==5 & gen.k>1) rn.mean1= abs.mean[yr.k,l.no.na,gen.k-1,scen.k,"rn"]
        
        if(gen.k==1) abs.mean1= abs.mean[yr.k,l.no.na,gen.k,scen.k,"absmid"]
        if(gen.k>1) abs.mean1= abs.mean[yr.k,l.no.na,gen.k-1,scen.k,"absmid"]
        
        #change NA values to negative values 
        abs.na.inds= abs.mean1[which( is.na(abs.mean1))]
        rn.na.inds= rn.mean1[which( is.na(rn.mean1))]
        
        #   #check abs mean
        #    if(!all(is.na(abs.mean1))){
        
        abs.mean1[which( is.na(abs.mean1))]= -10 
        rn.mean1[which( is.na(rn.mean1))]= -1000
        
        #Choose random sample of abs and rn values from current distribution (truncated normal) 
        abs.sample= sapply(abs.mean1, function(x) rtnorm(N.ind, mean = x, sd = abs.sd, lower=0.400, upper=0.700) )
        rn.sample= sapply(rn.mean1, function(x) rtnorm(N.ind, mean = x, sd = rn.sd, lower=-1, upper=1) )
        if(scen.mat[scen.k,1]==0) rn.sample[]=0
        
        #Add plasticity across sites and sample
        abs.plast <- abs.sample + rn.sample*(Tp-Tmid)
        #abs.mean[yr.k,,gen.k] <- abs.mean[yr.k,,gen.k]+abs.plast
        
        ##calculate fitness
        #use fitness function to predict Lambda for each individual
        #extract coefficients and calculate across abs samples
        fit.sample= foreach(cell.k=1:length(l.no.na), .combine="cbind") %do% {
          sapply(abs.plast[,cell.k], function(x) if( sum(is.na(fit[,cell.k]))==0) fit[1,cell.k]+x*fit[2,cell.k]+x^2*fit[3,cell.k] )
        } 
        #Fit.pred <- eval.fd(Abs.sample,Fitmod.year.gen) ### for spline
        
        #standardize to relative fitness and centered on trait mean
        fit.mean= colMeans(fit.sample)
        lambda.mean[yr.k,l.no.na,gen.k,scen.k]=fit.mean
        rel.fit= fit.sample/fit.mean
        
        absmid.dif= t( apply(abs.sample,1,'-',abs.mean1) )
        rn.dif= t( apply(rn.sample,1,'-',rn.mean1) )
        
        R2selnAbsmid<- rep(0, length(l.no.na) ) #No response to selection if no evolution
        R2selnRN<- rep(0, length(l.no.na) ) 
        #------------
        if(scen.k<5 & scen.mat[scen.k,2]==1){    
          ##selection analysis
          sel.fit= sapply(1:length(l.no.na), function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] +I(absmid.dif[,x]^2))$coefficients)
          
          #Save model
          sel.mod= sapply(1:length(l.no.na), function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] +I(absmid.dif[,x]^2) ) )
          ## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared
          
          #Response to selection
          BetaAbsmid <-sel.fit[2,]
          R2selnAbsmid <- h2*(abs.sd^2)*BetaAbsmid
        } #end scen.k<5
        #------------
        if(scen.k==5){    
          ##selection analysis
          sel.fit= sapply(1:length(l.no.na), function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] + rn.dif[,x] +I(absmid.dif[,x]^2) +I(rn.dif[,x]^2)+ rn.dif[,x]*absmid.dif[,x])$coefficients)
          
          #Save model
          sel.mod= sapply(1:length(l.no.na), function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] + rn.dif[,x] +I(absmid.dif[,x]^2) +I(rn.dif[,x]^2)+ rn.dif[,x]*absmid.dif[,x]) )
          ## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared
          
          #Response to selection
          BetaAbsmid <-sel.fit[2,]
          R2selnAbsmid <- h2*(abs.sd^2)*BetaAbsmid
          
          BetaRN <- sel.fit[3,] 
          R2selnRN <- h2*(rn.sd^2)*BetaRN
        } #end scen.k==5
        #-------------
        
        #Response to selection
        if(gen.k<3) {
          abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"absmid"]= abs.mean[yr.k,l.no.na,gen.k,scen.k,"absmid"] + R2selnAbsmid
          #Constain abs
          abs.mean[yr.k,which(abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"absmid"]>0.7),gen.k+1,scen.k,"absmid"]=0.7
          abs.mean[yr.k,which(abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"absmid"]<0.4),gen.k+1,scen.k,"absmid"]=0.4   
          
          #rn evolution
          if(scen.k==5){
            abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"rn"]= abs.mean[yr.k,l.no.na,gen.k,scen.k,"rn"] + R2selnRN
            #Constain abs
            abs.mean[yr.k,which(abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"rn"]>1),gen.k+1,scen.k,"rn"]= 1
            abs.mean[yr.k,which(abs.mean[yr.k,l.no.na,gen.k+1,scen.k,"rn"]< -1),gen.k+1,scen.k,"rn"]= -1
          }
        } #end evolutionary scenarios
        
        #Account for missing lambdas
        if(length(abs.na.inds)>0)  R2selnAbsmid[abs.na.inds]=NA    
        if(length(rn.na.inds)>0)   R2selnRN[rn.na.inds]=NA  
        
        #also put in next year's slot
        abs.mean[yr.k+1,l.no.na,1,scen.k,"absmid"]= abs.mean[yr.k,l.no.na,gen.k,scen.k,"absmid"] + R2selnAbsmid
        #Constain abs
        abs.mean[yr.k+1,which(abs.mean[yr.k+1,l.no.na,1,scen.k,"absmid"]>0.7),1,scen.k,"absmid"]=0.7
        abs.mean[yr.k+1,which(abs.mean[yr.k+1,l.no.na,1,scen.k,"absmid"]<0.4),1,scen.k,"absmid"]=0.4 
        
        if(scen.k==5) abs.mean[yr.k+1,l.no.na,1,scen.k,"rn"]= abs.mean[yr.k,l.no.na,gen.k,scen.k,"rn"] + R2selnRN
        #Constain abs
        abs.mean[yr.k+1,which(abs.mean[yr.k+1,l.no.na,gen.k,scen.k,"rn"]>1),1,scen.k,"rn"]= 1
        abs.mean[yr.k+1,which(abs.mean[yr.k+1,l.no.na,gen.k,scen.k,"rn"]< -1),1,scen.k,"rn"]= -1
        
        #Store other metrics
        abs.mean[yr.k,l.no.na,gen.k,scen.k,"abssample"]= colMeans(abs.plast)
        abs.mean[yr.k,l.no.na,gen.k,scen.k,"Babsmid"]= BetaAbsmid
        if(scen.k==5) abs.mean[yr.k,l.no.na,gen.k,scen.k,"Brn"]= BetaRN
        
      } #end scen loop
      
    } #Check lambda values exist
    
  } #end generation
  print(yr.k)
} #end year 

#=====================================
#Save output

setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/OUT/") #mac version

#saveRDS(abs.mean, "absmean.abs")
#saveRDS(lambda.mean, "lambdamean.abs")

abs.mean <- readRDS("absmean.abs")
lambda.mean <- readRDS("lambdamean.abs")

#================================

