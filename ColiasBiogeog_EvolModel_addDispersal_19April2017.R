#Evolutionary Model across biogeography
library(msm)
library(foreach)
library(reshape)

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
library(arrayhelpers)
library(gridExtra)
library(akima) #for interpolation
#-----------------------------
    
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

ngens=3 #adjust by elevation?

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
#previous version
#setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\3gen_rds")

#Lambda <- readRDS( paste("lambda1_",projs[proj.k],".rds", sep="") )
#pup.temps <- readRDS( paste("PupTemps_",projs[proj.k],".rds", sep="") )
##recent versions
Lambda <- readRDS( paste("lambda1_May12_",projs[proj.k],".rds", sep="") )
pup.temps <- readRDS( paste("PupTemps_May12_",projs[proj.k],".rds", sep="") )

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

#----------------------------
#Plot number generations by elevation

Lambda.yr.gen= Lambda[yr.k, , , gen.k, ] #150 841   7   3   4
# colSums(!is.na(Lambda[, , 1, 3, 1])) #CHECK NA counts

lgens= Lambda[, , 3, 2, 1]
lgens= cbind(pts.sel, t(lgens) )
lgens$lambda= lgens[,16]

ggplot(lgens, aes(x=elev, y=lambda)) +geom_point() 

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
    # colSums(!is.na(Lambda[, , 1, 3, 1])) #CHECK NA counts

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
    ## EXTRACT SUMMARY?:   fitr2 <- summary(lm.fitmod.yr)$r.squared
    
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
#_reg for regression 

setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")

#saveRDS(abs.mean, "absmean.abs")
#saveRDS(lambda.mean, "lambdamean.abs")

abs.mean <- readRDS("absmean.abs")
lambda.mean <- readRDS("lambdamean.abs")

#=====================================================
##  PLOT OUT

#sample sites to faciliate visualization
site.ind=sort(base::sample(1:nrow(pts.sel),200))

#=====================================================

#Fig 1. elev vs year for Jadult, Tpup, Tad    
  
  #dim(pup.temps)= 12 150 841 3

#Assess number generations
phen= cbind(pts.sel, t(pup.temps["Jadult",,,3]) ) 
phen= gather(phen, "year", "Jadult",9:(length(years)+8))
#count of Jadult=273
length(which(phen$Jadult==273))/nrow(phen)
#change 273 to NA
phen$Jadult[which(phen$Jadult==273)]=NA
#agregate sites that can't complete third generation
phen2= aggregate(phen, list(phen$ind), FUN=mean, na.rm=TRUE)
plot(phen2$elev, phen2$Jadult)

#Loop generations
for(gen in 1:3){

  #Jadult by gen
  phen= cbind(pts.sel, t(pup.temps["Jadult",inds,,gen]) ) 
  #subset points
  #phen= phen[site.ind,]
  
  phen= gather(phen, "year", "Jadult",9:(length(years)+8))
  phen$year= years[as.numeric(phen$year)]
  
  #replace Jadult=273 with NA
  phen$Jadult[which(phen$Jadult==273)]=NA
  #remove NAs
  phen= phen[which(!is.na(phen$Jadult)),]
  
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
  
  #plot single years and elevs to test
#  phen2= phen[phen$year==1960, ]
#  plot(phen2$elev, phen2$Jadult)
  
  #scatter plot
#  plot(years, phen[160,9:158])
#  plot(pts.sel$elev, phen[,130])
  
  #------------------
  
  #sample sites to faciliate visualization
  s.inds=sort(base::sample(1:nrow(phen),50000))  
  phen1= phen[s.inds,]
  
  #Interpolate
  s=interp(phen1$year,phen1$elev,phen1$Jadult, duplicate="strip")
  
  gdat <- interp2xyz(s, data.frame=TRUE)
  
  plot.Jad= ggplot(gdat) + 
    aes(x = x, y = y, z = z, fill = z) + 
    geom_tile() + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="Jadult") +
    theme_bw(base_size=16)+xlab("year")+ylab("elevation (m)")+theme(legend.position="bottom")
 
  #============================================================================
  #Tpup
  phen= cbind(pts.sel, t(pup.temps["Tpup",inds,,gen]) ) 
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
  
  #------------------
  
  #sample sites to faciliate visualization
  s.inds=sort(base::sample(1:nrow(phen),10000))  
  phen1= phen[s.inds,]
  
  #Interpolate
  s=interp(phen1$year,phen1$elev,phen1$Tpup, duplicate="strip")
  
  gdat <- interp2xyz(s, data.frame=TRUE)
  
  plot.Tpup= ggplot(gdat) + 
    aes(x = x, y = y, z = z, fill = z) + 
    geom_tile() + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="Tpup") +
    theme_bw(base_size=16)+xlab("year")+ylab("elevation (m)")+theme(legend.position="bottom")
  
  #------------------------------------
  #Tad
  phen= cbind(pts.sel, t(pup.temps["Tad",inds,,gen]) ) 
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
  
  #------------------
  
  #sample sites to faciliate visualization
  s.inds=sort(base::sample(1:nrow(phen),50000))  
  phen1= phen[s.inds,]
  #cut coldest temp to facilitate color scale
  phen1= phen1[which(phen1$Tad>7),]
  
  #Interpolate
  s=interp(phen1$year,phen1$elev,phen1$Tad, duplicate="strip")
  #, xo=seq(min(phen1$year), max(phen1$year), length = 80), yo=seq(min(phen1$elev), max(phen1$elev), length = 80))
  
  gdat <- interp2xyz(s, data.frame=TRUE)
  
  plot.Tad= ggplot(gdat) + 
    aes(x = x, y = y, z = z, fill = z) + 
    geom_tile() + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="Tad") +
    theme_bw(base_size=16)+xlab("year")+ylab("elevation (m)")+theme(legend.position="bottom")
  
  if(gen==1){plot.Jad1=plot.Jad; plot.Tpup1=plot.Tpup; plot.Tad1=plot.Tad}
  if(gen==2){plot.Jad2=plot.Jad; plot.Tpup2=plot.Tpup; plot.Tad2=plot.Tad}
  if(gen==3){plot.Jad3=plot.Jad; plot.Tpup3=plot.Tpup; plot.Tad3=plot.Tad}
   
} #end loop generations
  
  #------------------------------------
  # plot together
  setwd(paste(fdir,"figures\\", sep=""))
  
  pdf("Fig1_FigJadTpupTad.pdf", height=7, width=10)
  grid.arrange(plot.Jad1, plot.Tpup1, plot.Tad1,plot.Jad2, plot.Tpup2, plot.Tad2, ncol = 3, nrow=2)
  dev.off()
# plot.Jad3, plot.Tpup3, plot.Tad3,
  
#=======================================
  
#Fig 2. PLOT FITNESS CURVES
  #Lambda[years, sites, abs, gen, metrics: Lambda, FAT,Egg Viability]
  
  abs= seq(0.4,0.7,0.05)
  
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
  fc$elevation= cut(fc$elev, breaks=3)
  
  #restrict years
  fc= fc[which(fc$year %in% c(1980,2010,2040,2070) ),]
  fc$abs= abs[as.numeric(fc$abs)]  
  
  fc1= ddply(fc, .(elevation,year,abs), summarize, lambda=mean(lambda,na.rm=TRUE))
  fc1$year= as.factor(fc1$year)

  fcmap2 = ggplot(fc1, aes(x=abs, y=lambda, color=elevation, lty=year)) +geom_line(lwd=1) +theme_bw()+ theme(legend.position = "bottom")+ylim(1,2.4)
  
  #-------------------
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Fig2_FitnessCurves_elevs.pdf", height = 8, width = 8)
  
  print(fcmap2)
  
  dev.off()  

  #=======================================
  
  #Fig 3. OPTIMAL ABSORPTIVITY
  
  #AVERAGE ACROSS TIME PERIODS
  #initialize with optimum value yrs 1950-1980, across generations
  abs.opt.init1 <- colMeans(abs.opt[which(years %in% 1951:1980),,], na.rm = TRUE)
  colnames(abs.opt.init1)= c("gen1","gen2","gen3")
  abs.opt.init1= cbind(pts.sel,abs.opt.init1)
  abs.opt.init1$period=19511980
  
  abs.opt.init2 <- colMeans(abs.opt[which(years %in% 1981:2010),,], na.rm = TRUE)
  colnames(abs.opt.init2)= c("gen1","gen2","gen3")
  abs.opt.init2= cbind(pts.sel,abs.opt.init2)
  abs.opt.init2$period=19812010
  
  abs.opt.init3 <- colMeans(abs.opt[which(years %in% 2011:2040),,], na.rm = TRUE)
  colnames(abs.opt.init3)= c("gen1","gen2","gen3")
  abs.opt.init3= cbind(pts.sel,abs.opt.init3)
  abs.opt.init3$period=20112040
  
  abs.opt.init4 <- colMeans(abs.opt[which(years %in% 2041:2070),,], na.rm = TRUE)
  colnames(abs.opt.init4)= c("gen1","gen2","gen3")
  abs.opt.init4= cbind(pts.sel,abs.opt.init4)
  abs.opt.init4$period=20412070
  
  abs.opt.init= rbind(abs.opt.init1, abs.opt.init2, abs.opt.init3, abs.opt.init4 )
  
  #----------
  abso <- melt(abs.opt.init, value.name='abs',variable.name='gen',  measure.vars=c("gen1","gen2","gen3"))
  #fix names
  names(abso)[which(names(abso)=="variable")]="gen"
  names(abso)[which(names(abso)=="value")]="abs"
  
  abso$period=as.factor(abso$period)
  
  #restrict to first two generations
  abso= abso[which(abso$gen %in% c("gen1","gen2")),]
  
  #----------------------------
  #abs by gen
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Fig3_Abs_optbyElev.pdf", height = 6, width = 8)
  ggplot(abso, aes(x=elev, y=abs, color=period) ) +geom_point(shape=1)+ facet_grid(. ~ gen) +geom_smooth(method=lm, se=FALSE)+theme_bw()+ theme(legend.position = "bottom")
  dev.off()
  
  #=======================================
  
  #Fig 4. Fitness change by scenario
  
#  abs.mean #dims: yr.k, cell.k, gen.k, scen.k:no plast, plast, only plast, metrics: abssample, absmid, rn, Babsmid, Brn)
#  lambda.mean= array(NA, dim=c(length(years),nrow(pts.sel), 3, 5)) #dims: yr.k, cell.k, gen.k, scen.k
# scen.k: plast0evol0, plast1evol0, plast0evol1, plast1evol1, plast1evol1rnevol1

#calculate fitness differences by scenario  
lambda.diff= lambda.mean
lambda.diff[,,,2]= lambda.diff[,,,2]-lambda.diff[,,,1] 
lambda.diff[,,,3]= lambda.diff[,,,3]-lambda.diff[,,,1] 
lambda.diff[,,,4]= lambda.diff[,,,4]-lambda.diff[,,,1] 
lambda.diff[,,,5]= lambda.diff[,,,5]-lambda.diff[,,,1] 

#separate generations
lgen1= lambda.diff[,,1,]
lgen2= lambda.diff[,,2,]
lgen3= lambda.diff[,,3,]

#drop 1st scenario
lgen1= lgen1[,,2:5]
lgen2= lgen2[,,2:5]
lgen3= lgen3[,,2:5]

#--------------
#flatten array
l1= array2df(lgen1)
colnames(l1)=c("lambda","year","site","scen")
l2= array2df(lgen2)
colnames(l2)=c("lambda","year","site","scen")
l3= array2df(lgen3)
colnames(l3)=c("lambda","year","site","scen")

#add elevation
l1$elev= pts.sel$elev[l1$site]
l2$elev= pts.sel$elev[l2$site]
l3$elev= pts.sel$elev[l3$site]

#add year
l1$year= years[l1$year]
l2$year= years[l2$year]
l3$year= years[l3$year]

#add scenario
scens= c("plast0evol0", "plast1evol0", "plast0evol1", "plast1evol1", "plast1evol1rnevol1")
l1$scen= scens[l1$scen+1]
l2$scen= scens[l2$scen+1]
l3$scen= scens[l3$scen+1]

#PLOT
#Specify generation
lg=l1

#add time period
lg$period=NA
#lg$period[which(lg$year>1950 & lg$year<=1980)]= 19511980
lg$period[which(lg$year>1980 & lg$year<=2010)]= 19812010
lg$period[which(lg$year>2010 & lg$year<=2040)]= 20112040
lg$period[which(lg$year>2040 & lg$year<=2070)]= 20412070
lg$period= as.factor(lg$period)

#remove l1$period==NA
lg= lg[which(!is.na(lg$period)),]

#order scenarios
lg$scen= factor(lg$scen, levels=c("plast1evol0", "plast0evol1", "plast1evol1", "plast1evol1rnevol1") )

#group
lg.per= ddply(lg, .(period,elev,site,scen), summarize, lambda=mean(lambda,na.rm=TRUE))

#plot
lambda.scen= ggplot(lg.per)+aes(x = elev, y = lambda, color=period)+geom_point()+ facet_grid(. ~ scen)+ theme(legend.position = "bottom")+theme_bw()+ylim(-0.5,1)
#check ylim

#---------
#2nd gen
#Specify generation
lg=l2

#add time period
lg$period=NA
#lg$period[which(lg$year>1950 & lg$year<=1980)]= 19511980
lg$period[which(lg$year>1980 & lg$year<=2010)]= 19812010
lg$period[which(lg$year>2010 & lg$year<=2040)]= 20112040
lg$period[which(lg$year>2040 & lg$year<=2070)]= 20412070
lg$period= as.factor(lg$period)

#remove l1$period==NA
lg= lg[which(!is.na(lg$period)),]

#order scenarios
lg$scen= factor(lg$scen, levels=c("plast1evol0", "plast0evol1", "plast1evol1", "plast1evol1rnevol1") )

#group
lg.per= ddply(lg, .(period,elev,site,scen), summarize, lambda=mean(lambda,na.rm=TRUE))

#plot
lambda.scen2= ggplot(lg.per)+aes(x = elev, y = lambda, color=period)+geom_point()+ facet_grid(. ~ scen)+ theme(legend.position = "bottom")+theme_bw()+ylim(-0.5,1)

#====================================
#Compare absorptivities

#calculate absorptivity differences by scenario  
lambda.diff= abs.mean[,,,,2] #use absmid
lambda.diff[,,,2]= lambda.diff[,,,2]-lambda.diff[,,,1] 
lambda.diff[,,,3]= lambda.diff[,,,3]-lambda.diff[,,,1] 
lambda.diff[,,,4]= lambda.diff[,,,4]-lambda.diff[,,,1] 
lambda.diff[,,,5]= lambda.diff[,,,5]-lambda.diff[,,,1] 

#separate generations
lgen1= lambda.diff[,,1,]
lgen2= lambda.diff[,,2,]
lgen3= lambda.diff[,,3,]

#drop 1st scenario
lgen1= lgen1[,,2:5]
lgen2= lgen2[,,2:5]
lgen3= lgen3[,,2:5]

#--------------
#flatten array
l1= array2df(lgen1)
colnames(l1)=c("abs","year","site","scen")
l2= array2df(lgen2)
colnames(l2)=c("abs","year","site","scen")
l3= array2df(lgen3)
colnames(l3)=c("abs","year","site","scen")

#add elevation
l1$elev= pts.sel$elev[l1$site]
l2$elev= pts.sel$elev[l2$site]
l3$elev= pts.sel$elev[l3$site]

#add year
l1$year= years[l1$year]
l2$year= years[l2$year]
l3$year= years[l3$year]

#add scenario
scens= c("plast0evol0", "plast1evol0", "plast0evol1", "plast1evol1", "plast1evol1rnevol1")
l1$scen= scens[l1$scen+1]
l2$scen= scens[l2$scen+1]
l3$scen= scens[l3$scen+1]

#PLOT
#Specify generation
lg=l1

#add time period
lg$period=NA
#lg$period[which(lg$year>1950 & lg$year<=1980)]= 19511980
lg$period[which(lg$year>1980 & lg$year<=2010)]= 19812010
lg$period[which(lg$year>2010 & lg$year<=2040)]= 20112040
lg$period[which(lg$year>2040 & lg$year<=2070)]= 20412070
lg$period= as.factor(lg$period)

#order scenarios
lg$scen= factor(lg$scen, levels=c("plast1evol0", "plast0evol1", "plast1evol1", "plast1evol1rnevol1") )

#remove l1$period==NA
lg= lg[which(!is.na(lg$period)),]

#group
lg.per= ddply(lg, .(period,elev,site,scen), summarize, abs=mean(abs,na.rm=TRUE))

#plot
abs.scen= ggplot(lg.per)+aes(x = elev, y = abs, color=period)+geom_point()+ facet_grid(. ~ scen)+ theme(legend.position = "bottom")+theme_bw()+ylim(-0.2,0.2)
#check ylim

#-----------------------
#Lambda and abs by scenario

setwd(paste(fdir,"figures\\",sep="") )
pdf("Fig4_Lambda_Abs_scen.pdf", height = 12, width = 12)

grid.newpage()
pushViewport(viewport(layout=grid.layout(3,1)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

print(lambda.scen,vp=vplayout(1,1))
print(lambda.scen2,vp=vplayout(2,1))
print(abs.scen,vp=vplayout(3,1))

dev.off()

#==================================================
#PLOT ABS TIME SERIES

#Absorptivities across time and elevations

for(scen.k in 1:5){
  gen.k=1
  abs.all= cbind(pts.sel, t(abs.mean[,,gen.k,scen.k,"abssample"]) ) 
  abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
  abs.dat$year= years[as.numeric(abs.dat$year)]
  abs.dat$ecut= cut(abs.dat$elev, breaks=3)
  abs.agg1= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean, na.rm=TRUE)
  abs.dat1= abs.dat
  
  gen.k=2
  abs.all= cbind(pts.sel, t(abs.mean[,,gen.k,scen.k,"abssample"]) ) 
  abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
  abs.dat$year= years[as.numeric(abs.dat$year)]
  abs.dat$ecut= cut(abs.dat$elev, breaks=3)
  abs.agg2= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean, na.rm=TRUE)
  abs.dat2= abs.dat
  
  gen.k=3
  abs.all= cbind(pts.sel, t(abs.mean[,,gen.k,scen.k,"abssample"]) ) 
  abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
  abs.dat$year= years[as.numeric(abs.dat$year)]
  abs.dat$ecut= cut(abs.dat$elev, breaks=3)
  abs.agg3= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean, na.rm=TRUE)
  abs.dat3= abs.dat
  
  p.abs1all = ggplot(abs.dat1, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))#+ylim(0.6,0.75)
  p.abs2all = ggplot(abs.dat2, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))#+ylim(0.6,0.75)
  p.abs3all = ggplot(abs.dat3, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))#+ylim(0.6,0.75)
  
  p.abs1 = ggplot(abs.agg1, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))#+ylim(0.6,0.75)
  p.abs2 = ggplot(abs.agg2, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))#+ylim(0.6,0.75)
  p.abs3 = ggplot(abs.agg3, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))#+ylim(0.6,0.75)
  
  #Abs without plasticity
  gen.k=1
  abs.all= cbind(pts.sel, t(abs.mean[,,gen.k,scen.k,"absmid"]) ) 
  abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
  abs.dat$year= years[as.numeric(abs.dat$year)]
  abs.dat$ecut= cut(abs.dat$elev, breaks=3)
  abs.agg1m= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean, na.rm=TRUE)
  abs.dat1m= abs.dat
  
  gen.k=2
  abs.all= cbind(pts.sel, t(abs.mean[,,gen.k,scen.k,"absmid"]) ) 
  abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
  abs.dat$year= years[as.numeric(abs.dat$year)]
  abs.dat$ecut= cut(abs.dat$elev, breaks=3)
  abs.agg2m= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean, na.rm=TRUE)
  abs.dat2m= abs.dat
  
  gen.k=3
  abs.all= cbind(pts.sel, t(abs.mean[,,gen.k,scen.k,"absmid"]) ) 
  abs.dat= gather(abs.all, "year", "abs",9:(length(years)+8) )
  abs.dat$year= years[as.numeric(abs.dat$year)]
  abs.dat$ecut= cut(abs.dat$elev, breaks=3)
  abs.agg3m= aggregate(abs.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean, na.rm=TRUE)
  abs.dat3m= abs.dat
  
  p.abs1allm = ggplot(abs.dat1, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))#+ylim(0.5,0.80)
  p.abs2allm = ggplot(abs.dat2, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))#+ylim(0.5,0.80)
  p.abs3allm = ggplot(abs.dat3, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10))#+ylim(0.5,0.80)
  
  p.abs1m = ggplot(abs.agg1m, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10)) #+ylim(0.6,0.70)
  p.abs2m = ggplot(abs.agg2m, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10)) #+ylim(0.6,0.70)
  p.abs3m = ggplot(abs.agg3m, aes(x=year, y=abs, group=X, color=elev )) +geom_line() +theme_bw()+scale_color_gradientn(colours=matlab.like(10)) #+ylim(0.6,0.70)
  
  #Lambdas across time and elevations
  gen.k=1
  lambda.all= cbind(pts.sel, t(lambda.mean[,,gen.k,scen.k]) )
  lambda.dat= gather(lambda.all, "year", "lambda",9:(length(years)+8) )
  lambda.dat$year= years[as.numeric(lambda.dat$year)]
  lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
  lambda.agg1= aggregate(lambda.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean,na.rm=TRUE)
  lambda.dat1= lambda.dat
  
  gen.k=2
  lambda.all= cbind(pts.sel, t(lambda.mean[,,gen.k,scen.k]) )
  lambda.dat= gather(lambda.all, "year", "lambda",9:(length(years)+8) )
  lambda.dat$year= years[as.numeric(lambda.dat$year)]
  lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
  lambda.agg2= aggregate(lambda.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean,na.rm=TRUE)
  lambda.dat2= lambda.dat
  
  gen.k=3
  lambda.all= cbind(pts.sel, t(lambda.mean[,,gen.k,scen.k]) )
  lambda.dat= gather(lambda.all, "year", "lambda",9:(length(years)+8) )
  lambda.dat$year= years[as.numeric(lambda.dat$year)]
  lambda.dat$ecut= cut(lambda.dat$elev, breaks=3)
  lambda.agg3= aggregate(lambda.dat, list(abs.dat$year,abs.dat$ecut), FUN=mean,na.rm=TRUE)
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
#FIG 5: ABS WITH PLASTICITY

setwd(paste(fdir,"figures\\",sep="") )
pdf("FigSX_Abs_year.pdf", height = 12, width = 12)

grid.newpage()
pushViewport(viewport(layout=grid.layout(4,3)))
vplayout<-function(x,y)
  viewport(layout.pos.row=x,layout.pos.col=y)

#print(p.a11,vp=vplayout(1,1))
print(p.a21,vp=vplayout(1,1))
print(p.a31,vp=vplayout(2,1))
print(p.a41,vp=vplayout(3,1))
print(p.a51,vp=vplayout(4,1))

#print(p.a12,vp=vplayout(1,2))
print(p.a22,vp=vplayout(1,2))
print(p.a32,vp=vplayout(2,2))
print(p.a42,vp=vplayout(3,2))
print(p.a52,vp=vplayout(4,2))

#print(p.a13,vp=vplayout(1,3))
print(p.a23,vp=vplayout(1,3))
print(p.a33,vp=vplayout(2,3))
print(p.a43,vp=vplayout(3,3))
print(p.a53,vp=vplayout(4,3))

dev.off()

#===============================================
#MAPS LAMBDA AND ABS RELATIVE TO BASELINE (2010?)

#LAMBDA
#average lambdas across time period
#first generation
lambda.diff= lambda.mean[,,1,]

#change erroneous values
lambda.diff[which(lambda.diff< (-60))]=NA
lambda.diff[which(lambda.diff> 112)]=NA

per1= colMeans(lambda.diff[which(years %in% 1951:1980),,], na.rm = FALSE, dims = 1)
per2= colMeans(lambda.diff[which(years %in% 1981:2010),,], na.rm = FALSE, dims = 1)
per3= colMeans(lambda.diff[which(years %in% 2011:2040),,], na.rm = FALSE, dims = 1)
per4= colMeans(lambda.diff[which(years %in% 2041:2070),,], na.rm = FALSE, dims = 1)

#translate to difference from 1981-2010 for no plasticity or evolution
lper1s= per1 -per1[,1]
lper2s= per2 -per1[,1]
lper3s= per3 -per1[,1]
lper4s= per4 -per1[,1]

#----------
#second generation

lambda.diff= lambda.mean[,,2,]
per1= colMeans(lambda.diff[which(years %in% 1951:1980),,], na.rm = FALSE, dims = 1)
per2= colMeans(lambda.diff[which(years %in% 1981:2010),,], na.rm = FALSE, dims = 1)
per3= colMeans(lambda.diff[which(years %in% 2011:2040),,], na.rm = FALSE, dims = 1)
per4= colMeans(lambda.diff[which(years %in% 2041:2070),,], na.rm = FALSE, dims = 1)

#translate to difference from 1981-2010 for no plasticity or evolution
lper1s.gen2= per1 -per1[,1]
lper2s.gen2= per2 -per1[,1]
lper3s.gen2= per3 -per1[,1]
lper4s.gen2= per4 -per1[,1]
#----------

#ABS
#average abs across time period
#first generation
lambda.diff= abs.mean[,,1,,2]#abs mid

#change erroneous values
lambda.diff[which(lambda.diff>0.8 | lambda.diff<0.2)]=NA

per1= colMeans(lambda.diff[which(years %in% 1951:1980),,], na.rm = FALSE, dims = 1)
per2= colMeans(lambda.diff[which(years %in% 1981:2010),,], na.rm = FALSE, dims = 1)
per3= colMeans(lambda.diff[which(years %in% 2011:2040),,], na.rm = FALSE, dims = 1)
per4= colMeans(lambda.diff[which(years %in% 2041:2070),,], na.rm = FALSE, dims = 1)

#translate to difference from 1981-2010 for no plasticity or evolution
aper1s= per1 -per2[,1]
aper2s= per2 -per2[,1]
aper3s= per3 -per2[,1]
aper4s= per4 -per2[,1]

#determine breaks
l.breaks= rbind(lper1s, lper2s, lper3s, lper4s)
l.breaks=quantile(l.breaks, probs=seq(0,1,0.1), na.rm=TRUE)

a.breaks= rbind(aper1s, aper2s, aper3s, aper4s)
a.breaks=quantile(a.breaks, probs=seq(0,1,0.1), na.rm=TRUE)
 
l.lab= round(l.breaks, digits=2)
a.lab= round(a.breaks, digits=2)

#-------------------------------------------
#MAP

for(scen.k in 1:5){
  
  #LAMBDA
  #Set up data
  lper1= cbind(pts.sel, lper1s[,scen.k])
  lper2= cbind(pts.sel, lper2s[,scen.k])
  lper3= cbind(pts.sel, lper3s[,scen.k])
  lper4= cbind(pts.sel, lper4s[,scen.k])
  
  names(lper1)[9]="lambda"
  names(lper2)[9]="lambda"
  names(lper3)[9]="lambda"
  names(lper4)[9]="lambda"
              
  #omit NAs
  lper1= lper1[which(!is.na(lper1$lambda)),]
  lper2= lper2[which(!is.na(lper2$lambda)),]
  lper3= lper3[which(!is.na(lper3$lambda)),]
  lper4= lper4[which(!is.na(lper4$lambda)),]
  
  #set up map
  bbox <- ggmap::make_bbox(lon, lat, lper1, f = 0.1)
  map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')
  map1=ggmap(map_loc, margins=FALSE)
  
  lper1.map<- map1 + geom_raster(data=lper1, aes(fill = lambda), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=l.breaks, labels=l.lab )+ coord_fixed() + theme(legend.position="right")
  lper2.map<- map1 + geom_raster(data=lper2, aes(fill = lambda), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=l.breaks, labels=l.lab )+ coord_fixed() + theme(legend.position="right")
  lper3.map<- map1 + geom_raster(data=lper3, aes(fill = lambda), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=l.breaks, labels=l.lab)+ coord_fixed() + theme(legend.position="right")
  lper4.map<- map1 + geom_raster(data=lper4, aes(fill = lambda), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=l.breaks, labels=l.lab)+ coord_fixed() + theme(legend.position="right")
  
  #-------------
  
  #ABS
  #Set up data
  aper1= cbind(pts.sel, aper1s[,scen.k])
  aper2= cbind(pts.sel, aper2s[,scen.k])
  aper3= cbind(pts.sel, aper3s[,scen.k])
  aper4= cbind(pts.sel, aper4s[,scen.k])
  
  names(aper1)[9]="abs"
  names(aper2)[9]="abs"
  names(aper3)[9]="abs"
  names(aper4)[9]="abs"
  
  #omit NAs
  aper1= aper1[which(!is.na(aper1$abs)),]
  aper2= aper2[which(!is.na(aper2$abs)),]
  aper3= aper3[which(!is.na(aper3$abs)),]
  aper4= aper4[which(!is.na(aper4$abs)),]
  
  #set up map
  bbox <- ggmap::make_bbox(lon, lat, aper1, f = 0.1)
  map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')
  map1=ggmap(map_loc, margins=FALSE)
  
  aper1.map<- map1 + geom_raster(data=aper1, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=a.breaks, labels=a.lab)+ coord_fixed() + theme(legend.position="right")
  aper2.map<- map1 + geom_raster(data=aper2, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=a.breaks, labels=a.lab)+ coord_fixed() + theme(legend.position="right")
  aper3.map<- map1 + geom_raster(data=aper3, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=a.breaks, labels=a.lab)+ coord_fixed() + theme(legend.position="right")
  aper4.map<- map1 + geom_raster(data=aper4, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=a.breaks, labels=a.lab)+ coord_fixed() + theme(legend.position="right")
  
  #----------
  if(scen.k==1) {lper1.1=lper1.map; aper1.1=aper1.map;lper2.1=lper2.map; aper2.1=aper2.map; lper3.1=lper3.map; aper3.1=aper3.map; lper4.1=lper4.map; aper4.1=aper4.map}
  if(scen.k==2) {lper1.2=lper1.map; aper1.2=aper1.map;lper2.2=lper2.map; aper2.2=aper2.map; lper3.2=lper3.map; aper3.2=aper3.map; lper4.2=lper4.map; aper4.2=aper4.map}
  if(scen.k==3) {lper1.3=lper1.map; aper1.3=aper1.map;lper2.3=lper2.map; aper2.3=aper2.map; lper3.3=lper3.map; aper3.3=aper3.map; lper4.3=lper4.map; aper4.3=aper4.map}
  if(scen.k==4) {lper1.4=lper1.map; aper1.4=aper1.map;lper2.4=lper2.map; aper2.4=aper2.map; lper3.4=lper3.map; aper3.4=aper3.map; lper4.4=lper4.map; aper4.4=aper4.map}
  if(scen.k==5) {lper1.5=lper1.map; aper1.5=aper1.map;lper2.5=lper2.map; aper2.5=aper2.map; lper3.5=lper3.map; aper3.5=aper3.map; lper4.5=lper4.map; aper4.5=aper4.map}
  
} #end scen loop

#-------------
#CHECK MAP

#library(fields)
#quilt.plot(aper2$lon,aper2$lat, aper2$abs)

#Try interpolation
aper2.map<- map1 + geom_raster(data=aper2, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=a.breaks, labels=a.lab)+ coord_fixed() + theme(legend.position="right")

aper2.map<- map1 + geom_raster(data=aper2, interpolate=TRUE, aes(fill = abs), alpha=0.5)+ coord_cartesian()

#----------------------------

#FIG. 5: ABS MAP
setwd(paste(fdir,"figures\\",sep="") )
pdf("Fig5_Abs_map.pdf", height = 4, width = 8)

grid_arrange_shared_legend(aper1.3,aper2.3,aper3.3,aper4.3,aper1.4,aper2.4,aper3.4,aper4.4,aper1.5,aper2.5,aper3.5,aper4.5, ncol = 4, nrow = 3)

dev.off()

# aper1.1,aper2.1,aper3.1,aper4.1,aper1.2,aper2.2,aper3.2,aper4.2,

#grid.newpage()
#pushViewport(viewport(layout=grid.layout(3,4)))
#vplayout<-function(x,y)
#  viewport(layout.pos.row=x,layout.pos.col=y)

#print(aper1.3,vp=vplayout(1,1))
#print(aper2.3,vp=vplayout(2,1))
#print(aper3.3,vp=vplayout(3,1))
#print(aper4.3,vp=vplayout(4,1))

#print(aper1.4,vp=vplayout(1,2))
#print(aper2.4,vp=vplayout(2,2))
#print(aper3.4,vp=vplayout(3,2))
#print(aper4.4,vp=vplayout(4,2))

#print(aper1.5,vp=vplayout(1,3))
#print(aper2.5,vp=vplayout(2,3))
#print(aper3.5,vp=vplayout(3,3))
#print(aper4.5,vp=vplayout(4,3))

#-------------
#FIG. 6: LAMBDA MAP

pdf("Fig6_Lambda_map.pdf", height = 4, width = 8)

grid_arrange_shared_legend(lper1.3,lper2.3,lper3.3,lper4.3,lper1.4,lper2.4,lper3.4,lper4.4,lper1.5,lper2.5,lper3.5,lper4.5, ncol = 4, nrow = 3)

dev.off()

#lper1.1,lper2.1,lper3.1,lper4.1,lper1.2,lper2.2,lper3.2,lper4.2, 
#================================

#DISPERSAL
#Data from Watt et al. 1977, 1979

#Packages: fishmove, dispkernals

#Dispersal kernals across lattices: http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00117.x/full
#https://www.r-bloggers.com/continuous-dispersal-on-a-discrete-lattice/
#https://rdrr.io/rforge/ecomodtools/man/LatticeTransitionProbs.html

#Analyze number generations

#pts.sel # 841   8
#Lambda[yr.k, , , gen.k, ] #150 841   7   3   4
#pup.temps[yr.k, , , gen.k, ] # 12 150 841   3

lamb= cbind(pts.sel,t(Lambda[, , 1, 1, 1]) )

lamb[,"2015"]
lamb[,"2023"]
lamb[,"2024"]

#--------------
grid_arrange_shared_legend <- function(..., ncol = length(list(...)), nrow = 1, position = c("right")) {
  
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
  
  # return gtable invisibly
  invisible(combined)
  
}

#-------------------
#CHECK NAs in pupal temeprature etc

#Check previous versions
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\3gen_rds")
Lambda1 <- readRDS( paste("lambda1_",projs[proj.k],".rds", sep="") )
pup.temps1 <- readRDS( paste("PupTemps_",projs[proj.k],".rds", sep="") )
lper1= cbind(pts.sel, pup.temps1[9,10,,2])

#Set up data
lper1= cbind(pts.sel, pup.temps[9,10,,2]) #adult temps

names(lper1)[9]="temp"

#set up map
bbox <- ggmap::make_bbox(lon, lat, lper1, f = 0.1)
map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')
map1=ggmap(map_loc, margins=FALSE)

lper1.map<- map1 + geom_raster(data=lper1, aes(fill = temp), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(20))
#+ scale_fill_gradientn(colours=matlab.like(10), breaks=l.breaks, labels=l.lab )+ coord_fixed() + theme(legend.position="right")

lper1.map

#------------------------
  