#Evolutionary Model across biogeography
library(msm)
library(foreach)
library(reshape)
    
fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\"

#pick projection
proj.k=2
projs=c("bcc-csm","ccsm4","gfdl")
 
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
rn.sd= 0.0083

#-------------------------------
#EVO MODEL
N.ind=1000
a= seq(0.4,0.7,0.05)
#for finding a with max fitness
a.fit= as.data.frame(seq(0.4,0.7,0.01))
names(a.fit)="a"

#Read points
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")
pts.sel= read.csv( paste("COpoints.csv", sep="") ) #_",projs[proj.k],"
  
#Read lambdas and pupal temps
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")

Lambda <- readRDS( paste("lambda1_",projs[proj.k],".rds", sep="") )
pup.temps <- readRDS( paste("PupTemps_",projs[proj.k],".rds", sep="") )

#Find years with calculations
counts= rowSums(is.na(pup.temps[6,,,1]))
#inds=1:(which.min(counts==0) -1)
inds=1:150
years= years[inds]

#phenology trends
#plot(years,pup.temps[5,, 100, 1])
#plot(years,Lambda[, 1, 3, 1,1])

#==============================================
#Calculate optimal absorptivity

#Lambda[years, sites, abs, gen, metrics: Lambda, FAT,Egg Viability]
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
    ## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared
    
    #find maxima lambda
    abs.opt[yr.k,,gen.k]= as.vector(array(unlist(sapply(fit.mod, function(x) if(!is.null(x))a.fit$a[which.max(predict.lm(x, a.fit))] )), dim=c(1, nrow(pts.sel)) ) )
    
    } #end check data
  } #end gen loop
} #end year loop

#save optimal Alphas
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")

#saveRDS(abs.opt, paste("abs.opt_",projs[proj.k],".rds", sep=""))
abs.opt <- readRDS( paste("abs.opt_",projs[proj.k],".rds", sep="") )

#***************************************
#compute initial AbsMean 
int_elev = 0.4226; slope_elev = 0.06517
Tmid = 20; slope_plast = -0.0083  #if Tmid=22.5, -0.006667;

elev_km= pts.sel$elev/1000
#fill in missing elevations: http://www.inside-r.org/packages/cran/rgbif/docs/elevation

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
#-------------------------------
scen.mat= rbind(c(0,0,0),c(1,0,0),c(0,1,0),c(1,1,0),c(1,1,1) )
colnames(scen.mat)= c("plast","evol","evolRN"  )
 
for(yr.k in 1:length(years)) {
  
  ##loop through generations in each year
  for(gen.k in 1:ngens) {

    BetaAbsmid=NA
    
    Lambda.yr.gen= Lambda[yr.k, , , gen.k, ]

  if(all(!is.na(Lambda.yr.gen))){ #FIX TO DEAL WITH NAs
    #Extract temperatures
    Tp= pup.temps["Tpup",yr.k, , gen.k]
    
    #--------------------------
    #Fitness models
    #Estimate fitness functions across cells
    fit= array(unlist(apply(Lambda.yr.gen[,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2))$coefficients)), dim=c(3, nrow(pts.sel)) )
    #Save model
    fit.mod= apply(Lambda.yr.gen[,,1], 1, function(x) if(sum(is.na(x))==0) lm(x~a+I(a^2)) )
    ## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared
    
    #find maxima lambda
    abs.max= as.vector(array(unlist(sapply(fit.mod, function(x) if(!is.null(x))a.fit$a[which.max(predict.lm(x, a.fit))] )), dim=c(1, nrow(pts.sel)) ) )
    
    #-------------------------
    # LOOP PLASTICITY SCENARIOS
    for(scen.k in 1:5){ #plast0evol0, plast1evol0, plast0evol1, plast1evol1, plast1evol1rnevol1
    
    if(scen.mat[scen.k,1]==1) rn.mean1= rep(slope_plast, nrow(pts.sel) )
    if(scen.mat[scen.k,1]==0) rn.mean1= rep(0, nrow(pts.sel) )
    if(scen.k==5 & gen.k==1) rn.mean1= abs.mean[yr.k,,gen.k,scen.k,"rn"]
    if(scen.k==5 & gen.k>1) rn.mean1= abs.mean[yr.k,,gen.k-1,scen.k,"rn"]
    
    if(gen.k==1) abs.mean1= abs.mean[yr.k,,gen.k,scen.k,"absmid"]
    if(gen.k>1) abs.mean1= abs.mean[yr.k,,gen.k-1,scen.k,"absmid"]
    
    #change NA values to negative values 
    abs.na.inds= abs.mean1[which( is.na(abs.mean1))]
    rn.na.inds= rn.mean1[which( is.na(rn.mean1))]
    
    #check abs mean
    if(!all(is.na(abs.mean1))){
    
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
    fit.sample= foreach(cell.k=1:nrow(pts.sel), .combine="cbind") %do% {
      sapply(abs.plast[,cell.k], function(x) if( sum(is.na(fit[,cell.k]))==0) fit[1,cell.k]+x*fit[2,cell.k]+x^2*fit[3,cell.k] )
      } 
    #Fit.pred <- eval.fd(Abs.sample,Fitmod.year.gen) ### for spline

    #standardize to relative fitness and centered on trait mean
    fit.mean= colMeans(fit.sample)
    lambda.mean[yr.k,,gen.k,scen.k]=fit.mean
    rel.fit= fit.sample/fit.mean
    
    absmid.dif= t( apply(abs.sample,1,'-',abs.mean1) )
    rn.dif= t( apply(rn.sample,1,'-',rn.mean1) )

    R2selnAbsmid<- rep(0, nrow(pts.sel) ) #No response to selection if no evolution
    R2selnRN<- rep(0, nrow(pts.sel) ) 

    if(scen.k<5 & scen.mat[scen.k,2]==1){    
      ##selection analysis
      sel.fit= sapply(1:841, function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] +I(absmid.dif[,x]^2))$coefficients)
      
      #Save model
      sel.mod= sapply(1:841, function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] +I(absmid.dif[,x]^2) ) )
      ## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared
      
      #Response to selection
      BetaAbsmid <-sel.fit[2,]
      R2selnAbsmid <- h2*(abs.sd^2)*BetaAbsmid
    } #end scen.k<5

    if(scen.k==5){    
      ##selection analysis
      sel.fit= sapply(1:841, function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] + rn.dif[,x] +I(absmid.dif[,x]^2) +I(rn.dif[,x]^2)+ rn.dif[,x]*absmid.dif[,x])$coefficients)
      
      #Save model
      sel.mod= sapply(1:841, function(x) if(sum(is.na(x))==0) lm(rel.fit[,x]~absmid.dif[,x] + rn.dif[,x] +I(absmid.dif[,x]^2) +I(rn.dif[,x]^2)+ rn.dif[,x]*absmid.dif[,x]) )
      ## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared
    
    #Response to selection
    BetaAbsmid <-sel.fit[2,]
    R2selnAbsmid <- h2*(abs.sd^2)*BetaAbsmid
    
    BetaRN <- sel.fit[3,] 
    R2selnRN <- h2*(rn.sd^2)*BetaRN
    } #end scen.k==5

#    
    
#Response to selection
if(gen.k<3) {
  abs.mean[yr.k,,gen.k+1,scen.k,"absmid"]= abs.mean[yr.k,,gen.k,scen.k,"absmid"] + R2selnAbsmid
  #Constain abs
  abs.mean[yr.k,which(abs.mean[yr.k,,gen.k+1,scen.k,"absmid"]>0.7),gen.k+1,scen.k,"absmid"]=0.7
  abs.mean[yr.k,which(abs.mean[yr.k,,gen.k+1,scen.k,"absmid"]<0.4),gen.k+1,scen.k,"absmid"]=0.4   
  
  #rn evolution
  if(scen.k==5){
    abs.mean[yr.k,,gen.k+1,scen.k,"rn"]= abs.mean[yr.k,,gen.k,scen.k,"rn"] + R2selnRN
    #Constain abs
    abs.mean[yr.k,which(abs.mean[yr.k,,gen.k+1,scen.k,"rn"]>1),gen.k+1,scen.k,"rn"]= 1
    abs.mean[yr.k,which(abs.mean[yr.k,,gen.k+1,scen.k,"rn"]< -1),gen.k+1,scen.k,"rn"]= -1
  }
} 

} #check NA abs.mean1

#Accoutn for missing lambdas
if(length(abs.na.inds)>0)  R2selnAbsmid[abs.na.inds]=NA    
if(length(rn.na.inds)>0)   R2selnRN[rn.na.inds]=NA  
        
#also put in next year's slot
abs.mean[yr.k+1,,1,scen.k,"absmid"]= abs.mean[yr.k,,gen.k,scen.k,"absmid"] + R2selnAbsmid
#Constain abs
abs.mean[yr.k+1,which(abs.mean[yr.k+1,,1,scen.k,"absmid"]>0.7),1,scen.k,"absmid"]=0.7
abs.mean[yr.k+1,which(abs.mean[yr.k+1,,1,scen.k,"absmid"]<0.4),1,scen.k,"absmid"]=0.4 

if(scen.k==5) abs.mean[yr.k+1,,1,scen.k,"rn"]= abs.mean[yr.k,,gen.k,scen.k,"rn"] + R2selnRN
#Constain abs
abs.mean[yr.k+1,which(abs.mean[yr.k+1,,gen.k,scen.k,"rn"]>1),1,scen.k,"rn"]= 1
abs.mean[yr.k+1,which(abs.mean[yr.k+1,,gen.k,scen.k,"rn"]< -1),1,scen.k,"rn"]= -1

#Store other metrics
abs.mean[yr.k,,gen.k,scen.k,"abssample"]= colMeans(abs.plast)
abs.mean[yr.k,,gen.k,scen.k,"Babsmid"]= BetaAbsmid
if(scen.k==5) abs.mean[yr.k,,gen.k,scen.k,"Brn"]= BetaRN

} #end scen loop

 } #end check NA lambda

  } #end generation
  print(yr.k)
} #end year 

#=====================================
#Save output
#_reg for regression 

setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")

#saveRDS(abs.mean, "absmean.abs")
#saveRDS(lambda.mean, "lambdamean.abs")

abs.mean <- readRDS("absmean.abs")
lambda.mean <- readRDS("lambdamean.abs")

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

#sample sites to faciliate visualization
site.ind=sort(base::sample(1:nrow(pts.sel),200))

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
pdf(paste("FitnessCurves_20002075_",projs[proj.k],".pdf", sep=""), height = 12, width = 15)

#******************************** 

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
phen= gather(phen, "year", "Jadult",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]
p11 = ggplot(phen, aes(x=year, y=Jadult, group=ind, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(125,275)
#only subset of sites?

phen= cbind(pts.sel, t(pup.temps["Jadult",inds,,2]) )
phen= gather(phen, "year", "Jadult",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]
p12 = ggplot(phen, aes(x=year, y=Jadult, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(125,275)

phen= cbind(pts.sel, t(pup.temps["Jadult",inds,,3]) )
phen= gather(phen, "year", "Jadult",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]
p13 = ggplot(phen, aes(x=year, y=Jadult, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(125,275)

#Tpupal
phen= cbind(pts.sel, t(pup.temps["Tpup",inds,,1]) )
phen= gather(phen, "year", "Tpup",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]
p21 = ggplot(phen, aes(x=year, y=Tpup, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)

phen= cbind(pts.sel, t(pup.temps["Tpup",inds,,2]) )
phen= gather(phen, "year", "Tpup",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]
p22 = ggplot(phen, aes(x=year, y=Tpup, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)

phen= cbind(pts.sel, t(pup.temps["Tpup",inds,,3]) )
phen= gather(phen, "year", "Tpup",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]
p23 = ggplot(phen, aes(x=year, y=Tpup, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)

#Tadult
phen= cbind(pts.sel, t(pup.temps["Tad",inds,,1]) )
phen= gather(phen, "year", "Tad",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]
p31 = ggplot(phen, aes(x=year, y=Tad, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)

phen= cbind(pts.sel, t(pup.temps["Tad",inds,,2]) )
phen= gather(phen, "year", "Tad",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]
p32 = ggplot(phen, aes(x=year, y=Tad, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)

phen= cbind(pts.sel, t(pup.temps["Tad",inds,,3]) )
phen= gather(phen, "year", "Tad",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]
p33 = ggplot(phen, aes(x=year, y=Tad, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)
#---------------

setwd(paste(fdir,"figures\\",sep="") )
pdf("PhenTemp_year.pdf", height = 12, width = 12)

grid.newpage()
pushViewport(viewport(layout=grid.layout(3,3)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p11,vp=vplayout(1,1))
print(p12,vp=vplayout(1,2))
print(p13,vp=vplayout(1,3))
print(p21,vp=vplayout(2,1))
print(p22,vp=vplayout(2,2))
print(p23,vp=vplayout(2,3))
print(p31,vp=vplayout(3,1))
print(p32,vp=vplayout(3,2))
print(p33,vp=vplayout(3,3))

dev.off()
#==================================

#Absorptivities across time and elevations

for(scen.k in 1:5){
gen.k=1
abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"abssample"]) ) 
abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
abs.dat$year= years[as.numeric(abs.dat$year)]
abs.dat$ecut= cut(abs.dat$elev, breaks=3)
abs.agg1= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
abs.dat1= abs.dat

gen.k=2
abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"abssample"]) ) 
abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
abs.dat$year= years[as.numeric(abs.dat$year)]
abs.dat$ecut= cut(abs.dat$elev, breaks=3)
abs.agg2= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
abs.dat2= abs.dat

gen.k=3
abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"abssample"]) ) 
abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
abs.dat$year= years[as.numeric(abs.dat$year)]
abs.dat$ecut= cut(abs.dat$elev, breaks=3)
abs.agg3= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
abs.dat3= abs.dat

p.abs1all = ggplot(abs.dat1, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)
p.abs2all = ggplot(abs.dat2, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)
p.abs3all = ggplot(abs.dat3, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)

p.abs1 = ggplot(abs.agg1, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)
p.abs2 = ggplot(abs.agg2, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)
p.abs3 = ggplot(abs.agg3, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)

#Abs without plasticity
gen.k=1
abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"absmid"]) ) 
abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
abs.dat$year= years[as.numeric(abs.dat$year)]
abs.dat$ecut= cut(abs.dat$elev, breaks=3)
abs.agg1m= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
abs.dat1m= abs.dat

gen.k=2
abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"absmid"]) ) 
abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
abs.dat$year= years[as.numeric(abs.dat$year)]
abs.dat$ecut= cut(abs.dat$elev, breaks=3)
abs.agg2m= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
abs.dat2m= abs.dat

gen.k=3
abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"absmid"]) ) 
abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
abs.dat$year= years[as.numeric(abs.dat$year)]
abs.dat$ecut= cut(abs.dat$elev, breaks=3)
abs.agg3m= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
abs.dat3m= abs.dat

p.abs1allm = ggplot(abs.dat1, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)
p.abs2allm = ggplot(abs.dat2, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)
p.abs3allm = ggplot(abs.dat3, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)

p.abs1m = ggplot(abs.agg1m, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.6,0.70)
p.abs2m = ggplot(abs.agg2m, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.6,0.70)
p.abs3m = ggplot(abs.agg3m, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.6,0.70)

#Lambdas across time and elevations
gen.k=1
lambda.all= cbind(pts.sel, t(lambda.mean[inds,,gen.k,scen.k]) )
lambda.dat= gather(lambda.all, "year", "lambda",9:(length(years)+8) )
lambda.dat$year= years[as.numeric(lambda.dat$year)]
lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
lambda.agg1= aggregate(lambda.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
lambda.dat1= lambda.dat

gen.k=2
lambda.all= cbind(pts.sel, t(lambda.mean[inds,,gen.k,scen.k]) )
lambda.dat= gather(lambda.all, "year", "lambda",9:(length(years)+8) )
lambda.dat$year= years[as.numeric(lambda.dat$year)]
lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
lambda.agg2= aggregate(lambda.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
lambda.dat2= lambda.dat

gen.k=3
lambda.all= cbind(pts.sel, t(lambda.mean[inds,,gen.k,scen.k]) )
lambda.dat= gather(lambda.all, "year", "lambda",9:(length(years)+8) )
lambda.dat$year= years[as.numeric(lambda.dat$year)]
lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
lambda.agg3= aggregate(lambda.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
lambda.dat3= lambda.dat

p.lambda1all = ggplot(lambda.dat1, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
p.lambda2all = ggplot(lambda.dat2, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
p.lambda3all = ggplot(lambda.dat3, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))

p.lambda1 = ggplot(lambda.agg1, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
p.lambda2 = ggplot(lambda.agg2, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
p.lambda3 = ggplot(lambda.agg3, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))

if(scen.k==1) {p.a11=p.abs1; p.l11= p.lambda1; p.a12=p.abs2; p.l12= p.lambda2; p.a13=p.abs3; p.l13= p.lambda3; p.m11=p.abs1m; p.m12=p.abs2m; p.m13=p.abs3m}
if(scen.k==2) {p.a21=p.abs1; p.l21= p.lambda1; p.a22=p.abs2; p.l22= p.lambda2; p.a23=p.abs3; p.l23= p.lambda3; p.m21=p.abs1m; p.m22=p.abs2m; p.m23=p.abs3m}
if(scen.k==3) {p.a31=p.abs1; p.l31= p.lambda1; p.a32=p.abs2; p.l32= p.lambda2; p.a33=p.abs3; p.l33= p.lambda3; p.m31=p.abs1m; p.m32=p.abs2m; p.m33=p.abs3m}
if(scen.k==4) {p.a41=p.abs1; p.l41= p.lambda1; p.a42=p.abs2; p.l42= p.lambda2; p.a43=p.abs3; p.l43= p.lambda3; p.m41=p.abs1m; p.m42=p.abs2m; p.m43=p.abs3m}
if(scen.k==5) {p.a51=p.abs1; p.l51= p.lambda1; p.a52=p.abs2; p.l52= p.lambda2; p.a53=p.abs3; p.l53= p.lambda3; p.m51=p.abs1m; p.m52=p.abs2m; p.m53=p.abs3m}

if(scen.k==1) {p.a11all=p.abs1all; p.l11all= p.lambda1all; p.a12all=p.abs2all; p.l12all= p.lambda2all; p.a13all=p.abs3all; p.l13all= p.lambda3all; p.m11all=p.abs1allm; p.m12all=p.abs2allm; p.m13all=p.abs3allm}
if(scen.k==2) {p.a21all=p.abs1all; p.221all= p.lambda1all; p.a22all=p.abs2all; p.l22all= p.lambda2all; p.a23all=p.abs3all; p.l23all= p.lambda3all; p.m21all=p.abs1allm; p.m22all=p.abs2allm; p.m23all=p.abs3allm}
if(scen.k==3) {p.a31all=p.abs1all; p.331all= p.lambda1all; p.a32all=p.abs2all; p.l32all= p.lambda2all; p.a33all=p.abs3all; p.l33all= p.lambda3all; p.m31all=p.abs1allm; p.m32all=p.abs2allm; p.m33all=p.abs3allm}
if(scen.k==4) {p.a41all=p.abs1all; p.441all= p.lambda1all; p.a42all=p.abs2all; p.l42all= p.lambda2all; p.a43all=p.abs3all; p.l43all= p.lambda3all; p.m41all=p.abs1allm; p.m42all=p.abs2allm; p.m43all=p.abs3allm}
if(scen.k==5) {p.a51all=p.abs1all; p.551all= p.lambda1all; p.a52all=p.abs2all; p.l52all= p.lambda2all; p.a53all=p.abs3all; p.l53all= p.lambda3all; p.m51all=p.abs1allm; p.m52all=p.abs2allm; p.m53all=p.abs3allm}

} #end scen loop

#-----------------------
#ABS WITH PLASTICITY

setwd(paste(fdir,"figures\\",sep="") )
pdf("Abs_year.pdf", height = 12, width = 12)

grid.newpage()
pushViewport(viewport(layout=grid.layout(5,3)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p.a11,vp=vplayout(1,1))
print(p.a21,vp=vplayout(2,1))
print(p.a31,vp=vplayout(3,1))
print(p.a41,vp=vplayout(4,1))
print(p.a51,vp=vplayout(5,1))

print(p.a12,vp=vplayout(1,2))
print(p.a22,vp=vplayout(2,2))
print(p.a32,vp=vplayout(3,2))
print(p.a42,vp=vplayout(4,2))
print(p.a52,vp=vplayout(5,2))

print(p.a13,vp=vplayout(1,3))
print(p.a23,vp=vplayout(2,3))
print(p.a33,vp=vplayout(3,3))
print(p.a43,vp=vplayout(4,3))
print(p.a53,vp=vplayout(5,3))

dev.off()

#--------
#all elevations

setwd(paste(fdir,"figures\\",sep="") )
pdf("Abs_year_AllElev.pdf", height = 12, width = 12)

grid.newpage()
pushViewport(viewport(layout=grid.layout(5,3)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p.a11all,vp=vplayout(1,1))
print(p.a21all,vp=vplayout(2,1))
print(p.a31all,vp=vplayout(3,1))
print(p.a41all,vp=vplayout(4,1))
print(p.a51all,vp=vplayout(5,1))

print(p.a12all,vp=vplayout(1,2))
print(p.a22all,vp=vplayout(2,2))
print(p.a32all,vp=vplayout(3,2))
print(p.a42all,vp=vplayout(4,2))
print(p.a52all,vp=vplayout(5,2))

print(p.a13all,vp=vplayout(1,3))
print(p.a23all,vp=vplayout(2,3))
print(p.a33all,vp=vplayout(3,3))
print(p.a43all,vp=vplayout(4,3))
print(p.a53all,vp=vplayout(5,3))

dev.off()
#-----------------------
#ABS WITHOUT PLASTICITY

setwd(paste(fdir,"figures\\",sep="") )
pdf("Absmid_year.pdf", height = 12, width = 12)

grid.newpage()
pushViewport(viewport(layout=grid.layout(5,3)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p.m11,vp=vplayout(1,1))
print(p.m21,vp=vplayout(2,1))
print(p.m31,vp=vplayout(3,1))
print(p.m41,vp=vplayout(4,1))
print(p.m51,vp=vplayout(5,1))

print(p.m12,vp=vplayout(1,2))
print(p.m22,vp=vplayout(2,2))
print(p.m32,vp=vplayout(3,2))
print(p.m42,vp=vplayout(4,2))
print(p.m52,vp=vplayout(5,2))

print(p.m13,vp=vplayout(1,3))
print(p.m23,vp=vplayout(2,3))
print(p.m33,vp=vplayout(3,3))
print(p.m43,vp=vplayout(4,3))
print(p.m53,vp=vplayout(5,3))

dev.off()
#-----------
setwd(paste(fdir,"figures\\",sep="") )
pdf("Absmid_year_allElev.pdf", height = 12, width = 12)

grid.newpage()
pushViewport(viewport(layout=grid.layout(5,3)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p.m11all,vp=vplayout(1,1))
print(p.m21all,vp=vplayout(2,1))
print(p.m31all,vp=vplayout(3,1))
print(p.m41all,vp=vplayout(4,1))
print(p.m51all,vp=vplayout(5,1))

print(p.m12all,vp=vplayout(1,2))
print(p.m22all,vp=vplayout(2,2))
print(p.m32all,vp=vplayout(3,2))
print(p.m42all,vp=vplayout(4,2))
print(p.m52all,vp=vplayout(5,2))

print(p.m13all,vp=vplayout(1,3))
print(p.m23all,vp=vplayout(2,3))
print(p.m33all,vp=vplayout(3,3))
print(p.m43all,vp=vplayout(4,3))
print(p.m53all,vp=vplayout(5,3))

dev.off()

#-----------------------

setwd(paste(fdir,"figures\\",sep="") )
pdf("Lambda_year.pdf", height = 12, width = 12)

grid.newpage()
pushViewport(viewport(layout=grid.layout(5,3)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p.l11,vp=vplayout(1,1))
print(p.l21,vp=vplayout(2,1))
print(p.l31,vp=vplayout(3,1))
print(p.l41,vp=vplayout(4,1))
print(p.l51,vp=vplayout(5,1))

print(p.l12,vp=vplayout(1,2))
print(p.l22,vp=vplayout(2,2))
print(p.l32,vp=vplayout(3,2))
print(p.l42,vp=vplayout(4,2))
print(p.l52,vp=vplayout(5,2))

print(p.l13,vp=vplayout(1,3))
print(p.l23,vp=vplayout(2,3))
print(p.l33,vp=vplayout(3,3))
print(p.l43,vp=vplayout(4,3))
print(p.l53,vp=vplayout(5,3))

dev.off()

#----------
setwd(paste(fdir,"figures\\",sep="") )
pdf("Lambda_year_allElev.pdf", height = 12, width = 12)

grid.newpage()
pushViewport(viewport(layout=grid.layout(5,3)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(p.l11all,vp=vplayout(1,1))
print(p.l21all,vp=vplayout(2,1))
print(p.l31all,vp=vplayout(3,1))
print(p.l41all,vp=vplayout(4,1))
print(p.l51all,vp=vplayout(5,1))

print(p.l12all,vp=vplayout(1,2))
print(p.l22all,vp=vplayout(2,2))
print(p.l32all,vp=vplayout(3,2))
print(p.l42all,vp=vplayout(4,2))
print(p.l52all,vp=vplayout(5,2))

print(p.l13all,vp=vplayout(1,3))
print(p.l23all,vp=vplayout(2,3))
print(p.l33all,vp=vplayout(3,3))
print(p.l43all,vp=vplayout(4,3))
print(p.l53all,vp=vplayout(5,3))

dev.off()

#DO:
#CALCULATE ABS SLOPE BY DECADE?

#-----------------------------
#EVOLUTION OF RN

gen.k=1
scen.k=5
 
  #------------------------
  #FIG X. Selection gradients
  
  #specify generation and scenario
  gen.k=1
  scen.k=4
  
  #Selection on Absorptivities across time and elevations
  for(scen.k in 3:5){
    gen.k=1
    abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"Babsmid"]) ) 
    abs.dat= gather(abs.all, "year", "Babs",9:145)
    abs.dat$year= years[as.numeric(abs.dat$year)]
    abs.dat$ecut= cut(abs.dat$elev, breaks=3)
    abs.agg1= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
    abs.dat1= abs.dat
    
    gen.k=2
    abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"Babsmid"]) ) 
    abs.dat= gather(abs.all, "year", "Babs",9:145)
    abs.dat$year= years[as.numeric(abs.dat$year)]
    abs.dat$ecut= cut(abs.dat$elev, breaks=3)
    abs.agg2= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
    abs.dat2= abs.dat
    
    gen.k=3
    abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"Babsmid"]) ) 
    abs.dat= gather(abs.all, "year", "Babs",9:145)
    abs.dat$year= years[as.numeric(abs.dat$year)]
    abs.dat$ecut= cut(abs.dat$elev, breaks=3)
    abs.agg3= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
    abs.dat3= abs.dat
    
    # p.abs1 = ggplot(abs.dat1, aes(x=year, y=Babs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10)) #+facet_wrap(~ecut)
    # p.abs2 = ggplot(abs.dat2, aes(x=year, y=Babs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
    # p.abs3 = ggplot(abs.dat3, aes(x=year, y=Babs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
  
    #ave by elevation   
    p.abs1 = ggplot(abs.agg1, aes(x=year, y=Babs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
    p.abs2 = ggplot(abs.agg2, aes(x=year, y=Babs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
    p.abs3 = ggplot(abs.agg3, aes(x=year, y=Babs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
    
    if(scen.k==5){
      gen.k=1
      abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"Brn"]) ) 
      abs.dat= gather(abs.all, "year", "Brn",9:145)
      abs.dat$year= years[as.numeric(abs.dat$year)]
      abs.dat$ecut= cut(abs.dat$elev, breaks=3)
      abs.agg1= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
      abs.dat1= abs.dat
      
      gen.k=2
      abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"Brn"]) ) 
      abs.dat= gather(abs.all, "year", "Brn",9:145)
      abs.dat$year= years[as.numeric(abs.dat$year)]
      abs.dat$ecut= cut(abs.dat$elev, breaks=3)
      abs.agg2= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
      abs.dat2= abs.dat
      
      gen.k=3
      abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"Brn"]) ) 
      abs.dat= gather(abs.all, "year", "Brn",9:145)
      abs.dat$year= years[as.numeric(abs.dat$year)]
      abs.dat$ecut= cut(abs.dat$elev, breaks=3)
      abs.agg3= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
      abs.dat3= abs.dat
      
      # p.rn1 = ggplot(abs.dat1, aes(x=year, y=Brn, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
      # p.rn2 = ggplot(abs.dat2, aes(x=year, y=Brn, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
      # p.rn3 = ggplot(abs.dat3, aes(x=year, y=Brn, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
    
      #ave by elev
      p.rn1 = ggplot(abs.agg1, aes(x=year, y=Brn, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
      p.rn2 = ggplot(abs.agg2, aes(x=year, y=Brn, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
      p.rn3 = ggplot(abs.agg3, aes(x=year, y=Brn, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
      }
    
    if(scen.k==3) {p.a31=p.abs1; p.a32=p.abs2; p.a33=p.abs3}
    if(scen.k==4) {p.a41=p.abs1; p.a42=p.abs2; p.a43=p.abs3}
    if(scen.k==5) {p.a51=p.abs1; p.a52=p.abs2; p.a53=p.abs3}
    
  } #end scen loop
  
  #---------------
  
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Babs_year.pdf", height = 10, width = 12)
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(3,3)))
  vplayout<-function(x,y)
    viewport(layout.pos.row=x,layout.pos.col=y)
  
  print(p.a31,vp=vplayout(1,1))
  print(p.a41,vp=vplayout(2,1))
  print(p.a51,vp=vplayout(3,1))
  print(p.a32,vp=vplayout(1,2))
  print(p.a42,vp=vplayout(2,2))
  print(p.a52,vp=vplayout(3,2))
  print(p.a33,vp=vplayout(1,3))
  print(p.a43,vp=vplayout(2,3))
  print(p.a53,vp=vplayout(3,3))

  dev.off()
  
  #--------------
  
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Brn_year.pdf", height = 3, width = 12)
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1,3)))
  vplayout<-function(x,y)
    viewport(layout.pos.row=x,layout.pos.col=y)
  
  print(p.rn1,vp=vplayout(1,1))
  print(p.rn2,vp=vplayout(1,2))
  print(p.rn3,vp=vplayout(1,3))
  
  dev.off()
  
  #==============================
  #MAP
  
  for(scen.k in 1:5){
    
    #Set up data
    lambda.all= pts.sel
    lambda.all$lambda2000= lambda.mean[which(years=='2000'),,gen.k,scen.k]
    lambda.all$lambda2075= lambda.mean[which(years=='2075'),,gen.k,scen.k]
    lambda.all$abs2000= abs.mean[which(years=='2000'),,gen.k,scen.k,"absmid"]
    lambda.all$abs2075= abs.mean[which(years=='2075'),,gen.k,scen.k,"absmid"]
    
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
    
    if(scen.k==1) {a2000.1=abs2000map; a2075.1=abs2075map; l2000.1=lambda2000map; l2075.1=lambda2075map;}
    if(scen.k==2) {a2000.2=abs2000map; a2075.2=abs2075map; l2000.2=lambda2000map; l2075.2=lambda2075map;}
    if(scen.k==3) {a2000.3=abs2000map; a2075.3=abs2075map; l2000.3=lambda2000map; l2075.3=lambda2075map;}
    if(scen.k==4) {a2000.4=abs2000map; a2075.4=abs2075map; l2000.4=lambda2000map; l2075.4=lambda2075map;}
    if(scen.k==5) {a2000.5=abs2000map; a2075.5=abs2075map; l2000.5=lambda2000map; l2075.5=lambda2075map;}
    
  } #end scen loop
  
  #-------------------
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Abs_map.pdf", height = 15, width = 8)
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(5,2)))
  vplayout<-function(x,y)
    viewport(layout.pos.row=x,layout.pos.col=y)
  
  print(a2000.1,vp=vplayout(1,1))
  print(a2075.1,vp=vplayout(1,2))
  print(a2000.2,vp=vplayout(2,1))
  print(a2075.2,vp=vplayout(2,2))
  print(a2000.3,vp=vplayout(3,1))
  print(a2075.3,vp=vplayout(3,2))
  print(a2000.4,vp=vplayout(4,1))
  print(a2075.4,vp=vplayout(4,2))
  print(a2000.5,vp=vplayout(5,1))
  print(a2075.5,vp=vplayout(5,2))
  
  dev.off()
  
  #-------------
  pdf("Lambda_map.pdf", height = 15, width = 8)
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(5,2)))
  vplayout<-function(x,y)
    viewport(layout.pos.row=x,layout.pos.col=y)
  
  print(l2000.1,vp=vplayout(1,1))
  print(l2075.1,vp=vplayout(1,2))
  print(l2000.2,vp=vplayout(2,1))
  print(l2075.2,vp=vplayout(2,2))
  print(l2000.3,vp=vplayout(3,1))
  print(l2075.3,vp=vplayout(3,2))
  print(l2000.4,vp=vplayout(4,1))
  print(l2075.4,vp=vplayout(4,2))
  print(l2000.5,vp=vplayout(5,1))
  print(l2075.5,vp=vplayout(5,2))
  
  dev.off()
  
 
  #------------------------------  
#plot optimal across generation
   
  inds=1:120  #137
  
  gen.k=1
  lambda.all= cbind(pts.sel, t(abs.opt[inds,,gen.k]) )
#  lambda.all= lambda.all[site.ind,] #subsample to clarify plot
  lambda.dat= gather(lambda.all, "year", "abs",9:128)
  lambda.dat$year= years[as.numeric(lambda.dat$year)]
  lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
  lambda.agg1= aggregate(lambda.dat, list(lambda.dat$year,lambda.dat$ecut), FUN=mean)
  lambda.dat1= lambda.dat
  
  gen.k=2
  lambda.all= cbind(pts.sel, t(abs.opt[inds,,gen.k]) )
  #lambda.all= lambda.all[site.ind,]
  lambda.dat= gather(lambda.all, "year", "abs",9:128)
  lambda.dat$year= years[as.numeric(lambda.dat$year)]
  lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
  lambda.agg2= aggregate(lambda.dat, list(lambda.dat$year,lambda.dat$ecut), FUN=mean)
  lambda.dat2= lambda.dat
  
  gen.k=3
  lambda.all= cbind(pts.sel, t(abs.opt[inds,,gen.k]) )
 # lambda.all= lambda.all[site.ind,]
  lambda.dat= gather(lambda.all, "year", "abs",9:128)
  lambda.dat$year= years[as.numeric(lambda.dat$year)]
  lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
  lambda.agg3= aggregate(lambda.dat, list(lambda.dat$year,lambda.dat$ecut), FUN=mean)
  lambda.dat3= lambda.dat
  
  p.lambda1 = ggplot(lambda.agg1, aes(x=year, y=abs, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
  p.lambda2 = ggplot(lambda.agg2, aes(x=year, y=abs, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
  p.lambda3 = ggplot(lambda.agg3, aes(x=year, y=abs, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
  
  #-------------
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Abs_opt.pdf", height = 6, width = 12)
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1,3)))
  vplayout<-function(x,y)
    viewport(layout.pos.row=x,layout.pos.col=y)
  
  print(p.lambda1,vp=vplayout(1,1))
  print(p.lambda2,vp=vplayout(1,2))
  print(p.lambda3,vp=vplayout(1,3))
  
  dev.off()
  
  #----------------------------
  #PLOT ALL ELEVATIONs
  
  p.lambda1 = ggplot(lambda.dat1, aes(x=year, y=abs, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
  p.lambda2 = ggplot(lambda.dat2, aes(x=year, y=abs, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
  p.lambda3 = ggplot(lambda.dat3, aes(x=year, y=abs, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
   
  #-------------
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Abs_opt_allElevs.pdf", height = 6, width = 12)
  
  grid.newpage()
  pushViewport(viewport(layout=grid.layout(1,3)))
  vplayout<-function(x,y)
    viewport(layout.pos.row=x,layout.pos.col=y)
  
  print(p.lambda1,vp=vplayout(1,1))
  print(p.lambda2,vp=vplayout(1,2))
  print(p.lambda3,vp=vplayout(1,3))
  
  dev.off()
  
  #-----------------------------
  # PLOT OPTIMAL WITH AVERAGING ACROSS INITIAL TIME PERIOD
  
  #initialize with optimum value yrs 1950-1960, across generations
  abs.opt.init <- colMeans(abs.opt[1:30,,])
  colnames(abs.opt.init)= c("gen1","gen2","gen3")
  abs.opt.init= cbind(pts.sel,abs.opt.init)
  
  abso <- melt(abs.opt.init, value.name='abs',variable.name='gen',  measure.vars=c("gen1","gen2","gen3"))
  #fix names
  names(abso)[which(names(abso)=="variable")]="gen"
  names(abso)[which(names(abso)=="value")]="abs"
  
  #abs by gen
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Abs_optbyElev.pdf", height = 6, width = 12)
  ggplot(abso, aes(x=elev, y=abs, color=lat) ) +geom_point(shape=1)+ facet_grid(. ~ gen) +geom_smooth(method=lm, se=FALSE)
  dev.off()
  
  #opt abs across generations
  pdf("Abs_opt_AcrossGen.pdf", height = 6, width = 6)
  ggplot(abso, aes(x=gen, y=abs, color=elev, group=X) ) +geom_line()
  dev.off()
  
#=========================================================
#UPDATED PLOTS: FITNESS SURFACE

  library(gridExtra)
  library(akima) #for interpolations 
  library(tidyr) #for gather
  library(ggplot2)
  library(akima)
  library(fields) #for image plot

  #plot with shared legend
  grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("bottom", "right")) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <- ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x) x + theme(legend.position="none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(position,
                       "bottom" = arrangeGrob(do.call(arrangeGrob, gl),
                                              legend,
                                              ncol = 1,
                                              heights = unit.c(unit(1, "npc") - lheight, lheight)),
                       "right" = arrangeGrob(do.call(arrangeGrob, gl),
                                             legend,
                                             ncol = 2,
                                             widths = unit.c(unit(1, "npc") - lwidth, lwidth)))
    grid.newpage()
    grid.draw(combined)
    
  }
  
  #--------------------------------
#Fig 1. elev vs year for Jadult, Tpup, Tad

#dim(pup.temps)= 12 150 841 3
  
  #Jadult 1st gen
  phen= cbind(pts.sel, t(pup.temps["Jadult",inds,,1]) ) 
  phen= gather(phen, "year", "Jadult",9:(length(years)+8))
  phen$year= years[as.numeric(phen$year)]
  
  #sort by elevation
  ord= order(pts.sel$elev)
  z1= pup.temps["Jadult",inds,ord,1]
  elevs= pts.sel$elev[ord]
  dups= which(duplicated(elevs))
  elevs[dups]=elevs[dups]+0.2 #adjust duplicated elevations
  dups= which(duplicated(elevs))
  elevs[dups]=elevs[dups]+0.2 #adjust duplicated elevations
  dups= which(duplicated(elevs))
  elevs[dups]=elevs[dups]+0.2 #adjust duplicated elevations
  
  image.plot(x = years, y = elevs, z=z1)
  
  
 # test.spline <- Tps(data.frame(phen$year,phen$elev), phen$Jadult)
#  new.grid <- predictSurface(test.spline, nx = 200, ny = 200)
#  image(new.grid)
 
  
#------------------
#  z1= pup.temps["Jadult",inds,,1]
#  y01= seq(min(pts.sel$elev), max(pts.sel$elev), length.out=length(years) )
#  fld1<- bicubic(x = years, y = pts.sel$elev, z=z1, x0=years, y0=y01 )
#  gdat <- interp2xyz(fld1, data.frame=TRUE)
  
#  plot.Jad= ggplot(gdat) + 
#    aes(x = x, y = y, z = z, fill = z) + 
#    geom_tile() + 
#    geom_contour(color = "white", alpha = 0.5) + 
#    scale_fill_distiller(palette="Spectral", na.value="white", name="Jadult") +
#    theme_bw(base_size=18)+xlab("year")+ylab("elevation (m)")+theme(legend.position="bottom")


  #------------------------------------
  #Tpup
  phen= cbind(pts.sel, t(pup.temps["Tpup",inds,,1]) ) 
  phen= gather(phen, "year", "Tpup",9:(length(years)+8))
  phen$year= years[as.numeric(phen$year)]
  
  #sort by elevation
  ord= order(pts.sel$elev)
  z1= pup.temps["Tpup",inds,ord,1]
  elevs= pts.sel$elev[ord]
  dups= which(duplicated(elevs))
  elevs[dups]=elevs[dups]+0.2 #adjust duplicated elevations
  dups= which(duplicated(elevs))
  elevs[dups]=elevs[dups]+0.2 #adjust duplicated elevations
  dups= which(duplicated(elevs))
  elevs[dups]=elevs[dups]+0.2 #adjust duplicated elevations
  
  image.plot(x = years, y = elevs, z=z1)
  
  #------------------------------------
  #Tad
  phen= cbind(pts.sel, t(pup.temps["Tad",inds,,1]) ) 
  phen= gather(phen, "year", "Tad",9:(length(years)+8))
  phen$year= years[as.numeric(phen$year)]
  
  #sort by elevation
  ord= order(pts.sel$elev)
  z1= pup.temps["Tad",inds,ord,1]
  elevs= pts.sel$elev[ord]
  dups= which(duplicated(elevs))
  elevs[dups]=elevs[dups]+0.2 #adjust duplicated elevations
  dups= which(duplicated(elevs))
  elevs[dups]=elevs[dups]+0.2 #adjust duplicated elevations
  dups= which(duplicated(elevs))
  elevs[dups]=elevs[dups]+0.2 #adjust duplicated elevations
  
  image.plot(x = years, y = elevs, z=z1)
  
  #------------------------------------
# plot together
#    fig6= grid_arrange_shared_legend(f1,f2,f3,f4,f5, ncol = 3, nrow = 2)  
    
  
  
  
  