#Evolutionary Model across biogeography
library(msm)
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
rn.sd= 0.0083

inds=1:137
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
int_elev = 0.4226; slope_elev = 0.06517
Tmid = 20; slope_plast = -0.0083  #if Tmid=22.5, -0.006667;

elev_km= pts.sel$elev/1000
#fill in missing elevations: http://www.inside-r.org/packages/cran/rgbif/docs/elevation

abs.init <- int_elev+ slope_elev*elev_km

#initialize with optimum value yrs 1950-1960, across generations
abs.init2 <- rowMeans(colMeans(abs.opt[1:10,, ]))

plot(elev_km, abs.init, ylim=range(0.5, 0.7))
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

    Lambda.yr.gen= Lambda[yr.k, , , gen.k, ]

    #if(!is.na(Lambda.yr.gen)){ #FIX TO DEAL WITH NAs
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
    
    #Choose random sample of abs and rn values from current distribution (truncated normal) 
    abs.sample= sapply(abs.mean1, function(x) rtnorm(N.ind, mean = x, sd = abs.sd, lower=0.400, upper=0.700) )
    rn.sample= sapply(rn.mean1, function(x) rtnorm(N.ind, mean = x, sd = rn.sd, lower=-1, upper=1) )
    if(scen.mat[scen.k,1]==0) rn.sample[]=0
    
    #Add plasticity across sites and sample
    abs.plast <- abs.sample + rn.sample*(Tp-Tmid)
    #abs.mean[yr.k,,gen.k] <- abs.mean[yr.k,,gen.k]+abs.plast
            
    ##calculate fitness
    #Choose random sample of abs values from current distribution (truncated normal) 
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

    R2selnAbsmid<- rep(0, nrow(pts.sel) )
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

#Response to selection
if(gen.k<3) {
  abs.mean[yr.k,,gen.k+1,scen.k,"absmid"]= abs.mean[yr.k,,gen.k,scen.k,"absmid"] + R2selnAbsmid
  #Constain abs
  abs.mean[yr.k,which(abs.mean[yr.k,,gen.k+1,scen.k,"absmid"]>0.7),gen.k+1,scen.k,"absmid"]=0.7
  abs.mean[yr.k,which(abs.mean[yr.k,,gen.k+1,scen.k,"absmid"]<0.4),gen.k+1,scen.k,"absmid"]=0.4
  
  if(scen.k==5){
    abs.mean[yr.k,,gen.k+1,scen.k,"rn"]= abs.mean[yr.k,,gen.k,scen.k,"rn"] + R2selnRN
    #Constain abs
    abs.mean[yr.k,which(abs.mean[yr.k,,gen.k+1,scen.k,"rn"]>1),gen.k+1,scen.k,"rn"]= 1
    abs.mean[yr.k,which(abs.mean[yr.k,,gen.k+1,scen.k,"rn"]< -1),gen.k+1,scen.k,"rn"]= -1
  }
}

#also put in next year's slot
abs.mean[yr.k+1,,1,scen.k,"absmid"]= abs.mean[yr.k,,1,scen.k,"absmid"] + R2selnAbsmid
#Constain abs
abs.mean[yr.k+1,which(abs.mean[yr.k+1,,1,scen.k,"absmid"]>0.7),1,scen.k,"absmid"]=0.7
abs.mean[yr.k+1,which(abs.mean[yr.k+1,,1,scen.k,"absmid"]<0.4),1,scen.k,"absmid"]=0.4

if(scen.k==5){
  abs.mean[yr.k+1,,1,scen.k,"rn"]= abs.mean[yr.k,,1,scen.k,"rn"] + R2selnRN
  #Constain abs
  abs.mean[yr.k+1,which(abs.mean[yr.k+1,,1,scen.k,"rn"]> 1),1,scen.k,"rn"]=1
  abs.mean[yr.k+1,which(abs.mean[yr.k+1,,1,scen.k,"rn"]< -1),1,scen.k,"rn"]=-1
}

#Store other metrics
abs.mean[yr.k,,gen.k,scen.k,"abssample"]= colMeans(abs.plast)
abs.mean[yr.k,,gen.k,scen.k,"Babsmid"]= BetaAbsmid
if(scen.k==5) abs.mean[yr.k,,gen.k,scen.k,"Brn"]= BetaRN

} #end scen loop

 # } #end check NA lambda

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

#specify generation and scenario
gen.k=1
scen.k=4

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
p11 = ggplot(phen, aes(x=year, y=Jadult, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(125,275)

phen= cbind(pts.sel, t(pup.temps["Jadult",inds,,2]) )
phen= gather(phen, "year", "Jadult",9:145)
phen$year= years[as.numeric(phen$year)]
p12 = ggplot(phen, aes(x=year, y=Jadult, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(125,275)

phen= cbind(pts.sel, t(pup.temps["Jadult",inds,,3]) )
phen= gather(phen, "year", "Jadult",9:145)
phen$year= years[as.numeric(phen$year)]
p13 = ggplot(phen, aes(x=year, y=Jadult, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(125,275)

#Tpupal
phen= cbind(pts.sel, t(pup.temps["Tpup",inds,,1]) )
phen= gather(phen, "year", "Tpup",9:145)
phen$year= years[as.numeric(phen$year)]
p21 = ggplot(phen, aes(x=year, y=Tpup, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)

phen= cbind(pts.sel, t(pup.temps["Tpup",inds,,2]) )
phen= gather(phen, "year", "Tpup",9:145)
phen$year= years[as.numeric(phen$year)]
p22 = ggplot(phen, aes(x=year, y=Tpup, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)

phen= cbind(pts.sel, t(pup.temps["Tpup",inds,,3]) )
phen= gather(phen, "year", "Tpup",9:145)
phen$year= years[as.numeric(phen$year)]
p23 = ggplot(phen, aes(x=year, y=Tpup, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)

#Tadult
phen= cbind(pts.sel, t(pup.temps["Tad",inds,,1]) )
phen= gather(phen, "year", "Tad",9:145)
phen$year= years[as.numeric(phen$year)]
p31 = ggplot(phen, aes(x=year, y=Tad, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)

phen= cbind(pts.sel, t(pup.temps["Tad",inds,,2]) )
phen= gather(phen, "year", "Tad",9:145)
phen$year= years[as.numeric(phen$year)]
p32 = ggplot(phen, aes(x=year, y=Tad, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(5,25)

phen= cbind(pts.sel, t(pup.temps["Tad",inds,,3]) )
phen= gather(phen, "year", "Tad",9:145)
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
inds=1:137

for(scen.k in 1:5){
gen.k=1
abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"abssample"]) ) 
abs.dat= gather(abs.all, "year", "abs",9:145)
abs.dat$year= years[as.numeric(abs.dat$year)]
abs.dat$ecut= cut(abs.dat$elev, breaks=3)
abs.agg1= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
abs.dat1= abs.dat

gen.k=2
abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"abssample"]) ) 
abs.dat= gather(abs.all, "year", "abs",9:145)
abs.dat$year= years[as.numeric(abs.dat$year)]
abs.dat$ecut= cut(abs.dat$elev, breaks=3)
abs.agg2= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
abs.dat2= abs.dat

gen.k=3
abs.all= cbind(pts.sel, t(abs.mean[inds,,gen.k,scen.k,"abssample"]) ) 
abs.dat= gather(abs.all, "year", "abs",9:145)
abs.dat$year= years[as.numeric(abs.dat$year)]
abs.dat$ecut= cut(abs.dat$elev, breaks=3)
abs.agg3= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
abs.dat3= abs.dat

#p.abs1 = ggplot(abs.dat1, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)
#p.abs2 = ggplot(abs.dat2, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)
#p.abs3 = ggplot(abs.dat3, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)

p.abs1 = ggplot(abs.agg1, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)
p.abs2 = ggplot(abs.agg2, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)
p.abs3 = ggplot(abs.agg3, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.80)

#Lambdas across time and elevations
gen.k=1
lambda.all= cbind(pts.sel, t(lambda.mean[inds,,gen.k,scen.k]) )
lambda.dat= gather(lambda.all, "year", "lambda",9:145)
lambda.dat$year= years[as.numeric(lambda.dat$year)]
lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
lambda.agg1= aggregate(lambda.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
lambda.dat1= lambda.dat

gen.k=2
lambda.all= cbind(pts.sel, t(lambda.mean[inds,,gen.k,scen.k]) )
lambda.dat= gather(lambda.all, "year", "lambda",9:145)
lambda.dat$year= years[as.numeric(lambda.dat$year)]
lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
lambda.agg2= aggregate(lambda.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
lambda.dat2= lambda.dat

gen.k=3
lambda.all= cbind(pts.sel, t(lambda.mean[inds,,gen.k,scen.k]) )
lambda.dat= gather(lambda.all, "year", "lambda",9:145)
lambda.dat$year= years[as.numeric(lambda.dat$year)]
lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
lambda.agg3= aggregate(lambda.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean)
lambda.dat3= lambda.dat

#p.lambda1 = ggplot(lambda.dat1, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
#p.lambda2 = ggplot(lambda.dat2, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
#p.lambda3 = ggplot(lambda.dat3, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))

p.lambda1 = ggplot(lambda.agg1, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
p.lambda2 = ggplot(lambda.agg2, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))
p.lambda3 = ggplot(lambda.agg3, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))

if(scen.k==1) {p.a11=p.abs1; p.l11= p.lambda1; p.a12=p.abs2; p.l12= p.lambda2; p.a13=p.abs3; p.l13= p.lambda3}
if(scen.k==2) {p.a21=p.abs1; p.l21= p.lambda1; p.a22=p.abs2; p.l22= p.lambda2; p.a23=p.abs3; p.l23= p.lambda3}
if(scen.k==3) {p.a31=p.abs1; p.l31= p.lambda1; p.a32=p.abs2; p.l32= p.lambda2; p.a33=p.abs3; p.l33= p.lambda3}
if(scen.k==4) {p.a41=p.abs1; p.l41= p.lambda1; p.a42=p.abs2; p.l42= p.lambda2; p.a43=p.abs3; p.l43= p.lambda3}
if(scen.k==5) {p.a51=p.abs1; p.l51= p.lambda1; p.a52=p.abs2; p.l52= p.lambda2; p.a53=p.abs3; p.l53= p.lambda3}

} #end scen loop

#-----------------------
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

#DO:
#CALCULATE ABS SLOPE BY DECADE?
#ASSUME START AT OPTIMAL ABS. TOO COLD?

#-----------------------------
#EVOLUTION OF RN

gen.k=1
scen.k=5
inds=1:137
 
  #------------------------
  #FIG X. Selection gradients
  
  #specify generation and scenario
  gen.k=1
  scen.k=4
  
  #Selection on Absorptivities across time and elevations
  inds=1:137
  
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
  
  #==============================================
  #Plot optimal scenarios
  
  #Lambda[years, sites, abs, gen, metrics: Lambda, FAT,Egg Viability]
  abs.opt= array(NA, dim=c(length(years),nrow(pts.sel), 3))  
  
  for(yr.k in 1:length(years)) {
    
    ##loop through generations in each year
    for(gen.k in 1:ngens) {
      
      Lambda.yr.gen= Lambda[yr.k, , , gen.k, ]
      
      #if(!is.na(Lambda.yr.gen)){ #FIX TO DEAL WITH NAs
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
      abs.opt[yr.k,,gen.k]= as.vector(array(unlist(sapply(fit.mod, function(x) if(!is.null(x))a.fit$a[which.max(predict.lm(x, a.fit))] )), dim=c(1, nrow(pts.sel)) ) )
  
    } #end gen loop
  } #end year loop

  #save optimal Alphas
  setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")
  
  #saveRDS(abs.opt, "abs.opt")
  
  abs.opt <- readRDS("abs.opt")
  
  #------------------------------  
#plot optimal across generation
   
  inds=1:137
  
  gen.k=1
  lambda.all= cbind(pts.sel, t(abs.opt[inds,,gen.k]) )
  lambda.dat= gather(lambda.all, "year", "lambda",9:145)
  lambda.dat$year= years[as.numeric(lambda.dat$year)]
  lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
  lambda.agg1= aggregate(lambda.dat, list(lambda.dat$year,lambda.dat$ecut), FUN=mean)
  lambda.dat1= lambda.dat
  
  gen.k=2
  lambda.all= cbind(pts.sel, t(abs.opt[inds,,gen.k]) )
  lambda.dat= gather(lambda.all, "year", "lambda",9:145)
  lambda.dat$year= years[as.numeric(lambda.dat$year)]
  lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
  lambda.agg2= aggregate(lambda.dat, list(lambda.dat$year,lambda.dat$ecut), FUN=mean)
  lambda.dat2= lambda.dat
  
  gen.k=3
  lambda.all= cbind(pts.sel, t(abs.opt[inds,,gen.k]) )
  lambda.dat= gather(lambda.all, "year", "lambda",9:145)
  lambda.dat$year= years[as.numeric(lambda.dat$year)]
  lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
  lambda.agg3= aggregate(lambda.dat, list(lambda.dat$year,lambda.dat$ecut), FUN=mean)
  lambda.dat3= lambda.dat
  
  p.lambda1 = ggplot(lambda.agg1, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
  p.lambda2 = ggplot(lambda.agg2, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
  p.lambda3 = ggplot(lambda.agg3, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
  
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
  
  p.lambda1 = ggplot(lambda.dat1, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
  p.lambda2 = ggplot(lambda.dat2, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
  p.lambda3 = ggplot(lambda.dat3, aes(x=year, y=lambda, group=X, color=elev )) +geom_smooth(method=loess,se=FALSE) +theme_bw()+scale_color_gradientn(colours=matlab.like(10))+ylim(0.5,0.70)
   
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
  #### START SIMULATION AT OPTIMAL ALPHA