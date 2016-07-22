 #Evolutionary Model across biogeography

library(foreach)
    
fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\"
 
#EVOLUTIONARY MODEL
#For each year
#Fit fitness curves
#Simulate individuals: pick trait value
#Calculate fitness
#Evolution: update trait value
## Lambda[yr.k, cell.k, abs.k, gen.k, 1]

ngens=3
years= 1950:2099

#Store abs values
##abs= array(NA, dim=c(length(years),dim(Lambda)[2],3))
h2= 0.4
abs.mean=0.55; #UPDATE
abs.sd= 0.062
#-------------------------------
#EVO MODEL
N.ind=1000
a= seq(0.4,0.7,0.05)
#for finding a with max fitness
a.fit= as.data.frame(seq(0.4,0.7,0.01))
names(a.fit)="a"

#Read points
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")
pts.sel= read.csv("COpoints.csv")

#Read lambdas and pupal temps
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")

Lambda <- readRDS("lambda.rds")
pup.temps <- readRDS("PupTemps.rds")

#phenology trends
#plot(years,pup.temps[5,, 100, 1])
#plot(years,Lambda[, 1, 3, 1,1])

#***************************************
#compute initial AbsMean 
int_elev = 0.4282; slope_elev = 0.064
Tmid = 22.5; slope_plast =-0.006667;

elev_km= pts.sel$elev/1000
#fill in missing elevations: http://www.inside-r.org/packages/cran/rgbif/docs/elevation

abs.init <- int_elev+ slope_elev*elev_km

#-----------------------
#Save values
abs.mean= array(NA, dim=c(length(years),nrow(pts.sel), 3))  #dims: yr.k, cell.k, 3 (no plast, plast, only plast)
abs.mean[1,,1]= abs.init

lambda.mean= matrix(NA, nrow=length(years), ncol=nrow(pts.sel))  #dims: yr.k, cell.k, 3 (no plast, plast, only plast)
#-------------------------------

for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {

    Lambda.yr.gen= Lambda[yr.k, , , gen.k, ]

    #if(!is.na(Lambda.yr.gen)){ #FIX TO DEAL WITH NAs
    #Extract temperatures
    Tp= pup.temps["Tpup",yr.k, , gen.k]
     
    #Add plasticity
    abs.plast <- abs.mean[yr.k,,gen.k] + slope_plast*(Tp-Tmid)
    #abs.mean[yr.k,,gen.k] <- abs.mean[yr.k,,gen.k]+abs.plast
    
    #Choose random sample of abs values from current distribution (truncated normal) 
    abs.sample= sapply(abs.plast, function(x) rtnorm(N.ind, mean = x, sd = abs.sd, lower=0.400, upper=0.700) )
    
    #Estimate fitness functions across cells
    fit= array(unlist(apply(Lambda.yr.gen[,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2))$coefficients)), dim=c(3, length(abs.plast)))
    #Save model
    fit.mod= apply(Lambda.yr.gen[,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2)) )
    ## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared
   
    #find maxima lambda
    abs.max= as.vector(array(unlist(sapply(fit.mod, function(x) if(!is.null(x))a.fit$a[which.max(predict.lm(x, a.fit))] )), dim=c(1, length(abs.plast)) ) )
            
    ##calculate fitness
    #Choose random sample of abs values from current distribution (truncated normal) 
    #use fitness function to predict Lambda for each individual
    #extract coefficients and calculate across abs samples
    fit.sample= foreach(cell.k=1:nrow(pts.sel), .combine="cbind") %do% {
      sapply(abs.sample[,cell.k], function(x) if( sum(is.na(fit[,cell.k]))==0) fit[1,cell.k]+x*fit[2,cell.k]+x^2*fit[3,cell.k] )
      } 
    #Fit.pred <- eval.fd(Abs.sample,Fitmod.year.gen) ### for spline
    
#standardize to relative fitness and centered on trait mean
fit.mean= colMeans(fit.sample)
lambda.mean[yr.k,]=fit.mean
rel.fit= fit.sample/fit.mean

abs.dif= abs.sample- abs.plast

##selection analysis
sel.fit= array(unlist(sapply(1:841, function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~abs.dif[,x]+I(abs.dif[,x]^2))$coefficients)), dim=c(3, length(abs.plast)))
#Save model
sel.mod= array(unlist(sapply(1:841, function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~abs.dif[,x]+I(abs.dif[,x]^2)) )), dim=c(3, length(abs.plast)))
## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared

#Response to selection
R2seln= h2*(abs.sd^2)* sel.fit[2,]
#Response to selection
if(gen.k<3) {
  abs.mean[yr.k,,gen.k+1]= abs.mean[yr.k,,gen.k] + R2seln
  #Constain abs
  abs.mean[yr.k,which(abs.mean[yr.k,,gen.k+1]>0.7),gen.k+1]=0.7
  abs.mean[yr.k,which(abs.mean[yr.k,,gen.k+1]<0.4),gen.k+1]=0.4
}

#also put in next year's slot
abs.mean[yr.k+1,,1]= abs.mean[yr.k,,1] + R2seln
#Constain abs
abs.mean[yr.k+1,which(abs.mean[yr.k+1,,1]>0.7),1]=0.7
abs.mean[yr.k+1,which(abs.mean[yr.k+1,,1]<0.4),1]=0.4

 # } #end check NA lambda

  } #end generation
  print(yr.k)
} #end year 

#=====================================
#Save output

setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")

#saveRDS(abs.mean, "absmean.abs")
#saveRDS(lambda.mean, "lambdamean.abs")

#=====================================================
##  PLOT OUT
library(ggplot2)
library(maptools)
library(tidyr)
library(plyr)
#for mapping
library(ggmap)
library(maps)
library(mapdata)
library(colorRamps)     # for matlab.like(...)
library(grid)

#Absorptivities across time and elevations
inds=1:137

abs.all= cbind(pts.sel, t(abs.mean[inds,,1]) )
abs.dat= gather(abs.all, "year", "abs",9:145)
abs.dat$year= years[as.numeric(abs.dat$year)]

p.abs = ggplot(abs.dat, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))

#Lambdas across time and elevations
lambda.all= cbind(pts.sel, t(lambda.mean[inds,]) )
lambda.dat= gather(lambda.all, "year", "lambda",9:145)
lambda.dat$year= years[as.numeric(lambda.dat$year)]

p.lambda = ggplot(lambda.dat, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))

#---------------

setwd(paste(fdir,"figures\\",sep="") )
pdf("LambdaAbs_year.pdf", height = 5, width = 10)

grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p.abs,vp=vplayout(1,1))
print(p.lambda,vp=vplayout(1,2))

dev.off()

#DO:
#CALCULATE ABS SLOPE BY DECADE?
#ASSUME START AT OPTIMAL ABS. TOO COLD?

#--------------------- 
#MAP

#----------
#Set up data
lambda.all= pts.sel
lambda.all$lambda2000= lambda.mean[which(years=='2000'),]
lambda.all$lambda2075= lambda.mean[which(years=='2075'),]
lambda.all$abs2000= abs.mean[which(years=='2000'),,2]
lambda.all$abs2075= abs.mean[which(years=='2075'),,2]

#LAMBDAS
#MAP 2000
#subset to lambda>1
dat2= subset(lambda.all, lambda.all$lambda2000>1)

#set up map
bbox <- ggmap::make_bbox(lon, lat, dat2, f = 0.1)
map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')
map1=ggmap(map_loc, margins=FALSE)

lambda2000map<- map1 + geom_raster(data=dat2, aes(fill = lambda2000), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10))+ coord_fixed() + theme(legend.position="bottom")

#MAP 2075
#subset to lambda>1
dat2= subset(lambda.all, lambda.all$lambda2075>1)

lambda2075map<- map1 + geom_raster(data=dat2, aes(fill = lambda2075), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10))+ coord_fixed() + theme(legend.position="bottom")
#----------
#ABSORPTIVITIES

dat2=lambda.all

#MAP 2000
abs2000map<- map1 + geom_raster(data=dat2, aes(fill = abs2000), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10))+ coord_fixed() + theme(legend.position="bottom")

#MAP 2075
abs2075map<- map1 + geom_raster(data=dat2, aes(fill = abs2075), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10))+ coord_fixed() + theme(legend.position="bottom")

#-------------------
setwd(paste(fdir,"figures\\",sep="") )
pdf("LambdaAbs_map.pdf", height = 12, width = 15)

grid.newpage()
pushViewport(viewport(layout=grid.layout(2,2)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(abs2000map,vp=vplayout(1,1))
print(abs2075map,vp=vplayout(1,2))
print(lambda2000map,vp=vplayout(2,1))
print(lambda2075map,vp=vplayout(2,2))

dev.off()

#=====================================================
#Fig X. PLOT FITNESS CURVES
#Lambda[years, sites, abs, gen, metrics: Lambda, FAT,Egg Viability]

abs= seq(0.4,0.7,0.05)

dat= Lambda[years %in% c(2000,2075),,,2,1]
dat= aperm(dat, perm=c(2,1,3))

dat1= cbind(pts.sel, dat[,1,])
dat2= cbind(pts.sel, dat[,2,])

fc2000= gather(dat1, "abs", "lambda",9:15)
fc2075= gather(dat2, "abs", "lambda",9:15)

fc2000$year="2000"
fc2075$year="2075"
fc2000$group=paste(fc2000$X,"2000",sep="")
fc2075$group=paste(fc2075$X,"2075",sep="")

fc=rbind(fc2000,fc2075)
fc$abs= abs[as.numeric(fc$abs)]

fcmap = ggplot(fc, aes(x=abs, y=lambda, group=group, color=elev)) +facet_wrap(~year) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
#-------------------
setwd(paste(fdir,"figures\\",sep="") )
pdf("FitnessCurves_20002075.pdf", height = 12, width = 15)

print(fcmap)

dev.off()

#=================================
#group by elevation
dat= Lambda[,,,2,1]

dat1= cbind(pts.sel, t(dat[,,1]) )
dat2= cbind(pts.sel, t(dat[,,2]) )
dat3= cbind(pts.sel, t(dat[,,3]) )
dat4= cbind(pts.sel, t(dat[,,4]) )
dat5= cbind(pts.sel, t(dat[,,5]) )
dat6= cbind(pts.sel, t(dat[,,6]) )
dat7= cbind(pts.sel, t(dat[,,7]) )

fc1= gather(dat1, "year", "lambda",9:158)
fc1$abs=1
fc2= gather(dat2, "year", "lambda",9:158)
fc2$abs=2
fc3= gather(dat3, "year", "lambda",9:158)
fc3$abs=3
fc4= gather(dat4, "year", "lambda",9:158)
fc4$abs=4
fc5= gather(dat5, "year", "lambda",9:158)
fc5$abs=5
fc6= gather(dat6, "year", "lambda",9:158)
fc6$abs=6
fc7= gather(dat7, "year", "lambda",9:158)
fc7$abs=7

fc=rbind(fc1,fc2,fc3,fc4,fc5,fc6,fc7)
fc$ecut= cut(fc$elev, breaks=3)

fc1= ddply(fc, .(ecut,year,abs), summarize, lambda=mean(lambda,na.rm=TRUE))
fc1$year= as.numeric(fc1$year)

fc1$abs= abs[as.numeric(fc1$abs)]     

fcmap2 = ggplot(fc1, aes(x=abs, y=lambda, group= as.factor(year),color=year)) +facet_wrap(~ecut) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))

#-------------------
setwd(paste(fdir,"figures\\",sep="") )
pdf("FitnessCurves_elevs.pdf", height = 12, width = 15)

print(fcmap2)

dev.off()

#==========================================
#FIG X. PHENOLOGY / TEMPS PLOTS

#pup.temps[stats, years, sites, gen] 
#STATS: pup.temps[[1]]= c("stat","yr","gen","Jlarv", "Jpup","Jadult","Tlarv","Tpup","Tad","Tlarv_fixed","Tpup_fixed","Tad_fixed") 

#Jadult
phen= cbind(pts.sel, t(pup.temps["Jadult",inds,,1]) )
phen= gather(phen, "year", "Jadult",9:145)
phen$year= years[as.numeric(phen$year)]

p1 = ggplot(phen, aes(x=year, y=Jadult, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))

#Tpupal
phen= cbind(pts.sel, t(pup.temps["Tpup",inds,,1]) )
phen= gather(phen, "year", "Tpup",9:145)
phen$year= years[as.numeric(phen$year)]

p2 = ggplot(phen, aes(x=year, y=Tpup, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))

#Tadult
phen= cbind(pts.sel, t(pup.temps["Tad",inds,,1]) )
phen= gather(phen, "year", "Tad",9:145)
phen$year= years[as.numeric(phen$year)]

p3 = ggplot(phen, aes(x=year, y=Tad, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
#---------------

setwd(paste(fdir,"figures\\",sep="") )
pdf("PhenTemp_year.pdf", height = 5, width = 15)

grid.newpage()
pushViewport(viewport(layout=grid.layout(1,3)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p1,vp=vplayout(1,1))
print(p2,vp=vplayout(1,2))
print(p3,vp=vplayout(1,3))

dev.off()


