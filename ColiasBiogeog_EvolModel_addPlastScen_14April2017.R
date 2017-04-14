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

#Fig 1. elev vs year for Jadult, Tpup, Tad    
  
  #dim(pup.temps)= 12 150 841 3
  
  #Jadult 1st gen
  phen= cbind(pts.sel, t(pup.temps["Jadult",inds,,1]) ) 
  #subset points
  #phen= phen[site.ind,]
  
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
  
  #plot single years and elevs to test
  phen2= phen[phen$year==1960, ]
  plot(phen2$elev, phen2$Jadult)
  
  #scatter plot
  plot(years, phen[160,9:158])
  plot(pts.sel$elev, phen[,130])
  
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
  #geom_contour(color = "white", alpha = 0.5) + 
  plot.Jad
  #============================================================================
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
  #geom_contour(color = "white", alpha = 0.5) + 
  plot.Tpup
  
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
  #geom_contour(color = "white", alpha = 0.5) + 
  plot.Tad
  
  #------------------------------------
  # plot together
  
  #fig6= grid_arrange_shared_legend(f1,f2,f3,f4,f5, ncol = 3, nrow = 2)  
  setwd(paste(fdir,"figures\\", sep=""))
  
  pdf("FigJadTpupTad.pdf", height=4, width=10)
  grid.arrange(plot.Jad, plot.Tpup, plot.Tad, ncol = 3)
  dev.off()

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
  fc= fc[which(fc$year %in% c(1980,2010,2040) ),]
  fc$abs= abs[as.numeric(fc$abs)]  
  
  fc1= ddply(fc, .(elevation,year,abs), summarize, lambda=mean(lambda,na.rm=TRUE))
  fc1$year= as.factor(fc1$year)

  fcmap2 = ggplot(fc1, aes(x=abs, y=lambda, color=elevation, lty=year)) +geom_line(lwd=1) +theme_bw()+ theme(legend.position = "bottom")
  
  #-------------------
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("FitnessCurves_elevs.pdf", height = 8, width = 8)
  
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
  
  abs.opt.init= rbind(abs.opt.init1, abs.opt.init2, abs.opt.init3 )
  
  #----------
  abso <- melt(abs.opt.init, value.name='abs',variable.name='gen',  measure.vars=c("gen1","gen2","gen3"))
  #fix names
  names(abso)[which(names(abso)=="variable")]="gen"
  names(abso)[which(names(abso)=="value")]="abs"
  
  abso$period=as.factor(abso$period)
  
  #----------------------------
  #abs by gen
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Fig3_Abs_optbyElev.pdf", height = 6, width = 12)
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

#group across years: 1980-2010, 2010-2040

#AVERAGE ACROSS TIME PERIODS
#initialize with optimum value yrs 1950-1980, across generations
lambda.diff2 <- colMeans(lambda.diff[which(years %in% 1981:2010),,,], na.rm = TRUE)
#generation 1
lambda.diff2 <- lambda.diff2[,1,]
lambda.diff2= cbind(pts.sel,lambda.diff2)
lambda.diff2$period= 19812010

lambda.diff3 <- colMeans(lambda.diff[which(years %in% 2011:2040),,,], na.rm = TRUE)
#generation 1
lambda.diff3 <- lambda.diff3[,1,]
lambda.diff3= cbind(pts.sel,lambda.diff3)
lambda.diff3$period= 20112040

lambda.diffs= rbind(lambda.diff2, lambda.diff3 )

#----------
ldiff <- melt(lambda.diffs, value.name='lambda',variable.name='scen',  measure.vars=c("1","2","3","4","5"))
#fix names
names(ldiff)[which(names(ldiff)=="variable")]="scen"
names(ldiff)[which(names(ldiff)=="value")]="lambda"

#cut scenerio 1
ldiff= ldiff[ldiff$scen!=1,]

ldiff$period= as.factor(ldiff$period)

scens= c("plast0evol0", "plast1evol0", "plast0evol1", "plast1evol1", "plast1evol1rnevol1")
ldiff$scen= scens[ldiff$scen]

#point plots
ggplot(ldiff)+aes(x = elev, y = lambda, color=period)+geom_point()+ facet_grid(. ~ scen)+ theme(legend.position = "bottom")

#--------------
#all years
lgen1= lambda.diff[,,1,]
lgen2= lambda.diff[,,2,]
lgen3= lambda.diff[,,3,]


l1= adply(lambda.diff, c(4)) 

l1= as.data.frame(as.table(lambda.diff))

l1= lambda.diff[,,1,]


ldiff <- melt(lambda.diff, value.name='lambda',variable.name='scen',  measure.vars=c("1","2","3","4","5"))


ggplot(lambda.diff)+aes(x = elev, y = lambda, color=period)+geom_point()+ facet_grid(. ~ scen)+ theme(legend.position = "bottom")





  
  