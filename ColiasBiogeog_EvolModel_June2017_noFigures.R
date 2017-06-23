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
library(scales)
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

#RESOURCES FOR ADDING DISPERSAL KERNAL
#USE:
#***  https://rdrr.io/rforge/ecomodtools/man/LatticeTransitionProbs.html

# https://www.r-bloggers.com/continuous-dispersal-on-a-discrete-lattice/
# http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00117.x/full
#useful? https://github.com/pedroj/dispkernels

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

Lambda <- readRDS( paste("lambda1_",projs[proj.k],".rds", sep="") )
pup.temps <- readRDS( paste("PupTemps_",projs[proj.k],".rds", sep="") )
##recent versions
#Lambda <- readRDS( paste("lambda1_May12_",projs[proj.k],".rds", sep="") )
#pup.temps <- readRDS( paste("PupTemps_May12_",projs[proj.k],".rds", sep="") )

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
#gen.k=3
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

#================================

#DISPERSAL
#Data from Watt et al. 1977, 1979

#Packages: fishmove, dispkernals

#Dispersal kernals across lattices: http://onlinelibrary.wiley.com/doi/10.1111/j.2041-210X.2011.00117.x/full
#https://www.r-bloggers.com/continuous-dispersal-on-a-discrete-lattice/

#RESOURCES FOR ADDING DISPERSAL KERNAL
#USE:
#***  https://rdrr.io/rforge/ecomodtools/man/LatticeTransitionProbs.html

#------------------------------------
#install.packages("ecomodtools", repos="http://R-Forge.R-project.org")
#library(ecomodtools)

# 2D Gaussian dispersal function defined according to Clark et al. (1999)
# Used in the Monte Carlo integration
genexp.disp <- function(from, to, params)
{
  r.sq <- (from[1] - to[1])^2 + (from[2] - to[2])^2
  (1.0 / (pi * params[1]^2)) * exp(-(r.sq / params[1]^2))
}

# Monte Carlo integration settings
max.prob <- .Machine$double.eps^0.25
initial.step <-1 #10000
step.size <-1 #10000
max.runs <-100 #1000000
max.rel.var <- .Machine$double.eps^0.25

# Boundary condition
boundary <- "restricting"

# Alpha parameter setting
params <- c(1.0)

# Make a small grid (3 x 3) for testing purposes
x.coords <- seq(0.0, 3.0, 1.0)
y.coords <- seq(0.0, 3.0, 1.0)
test.grid <- MakeGridFromCoords(x.coords, y.coords)

##Manually make grid
x1= c(1,2,3,1,2,3,1,2,3)
x2= c(2,3,4,2,3,4,2,3,4)
y1= c(1,1,1,2,2,2,3,3,3)
y2= c(2,2,2,3,3,3,4,4,4)
test.grid= as.data.frame(cbind(x1,x2,y1,y2))


# Calculate the transition matrices for each different approximation method

# Centroid-to-centroid transition matrix using Monte Carlo integration
CC.mc <- LatticeTransitionProbs(
  x1 = test.grid$x1, x2 = test.grid$x2, y1 = test.grid$y1, y2 = test.grid$y2,
  func = genexp.disp, approx.method = "CC", boundary = boundary, max.prob = max.prob,
  initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs, params = params)
# Centroid-to-centroid transition matrix using analytic results
CC.an <- LatticeTransitionProbs(
  x1 = test.grid$x1, x2 = test.grid$x2, y1 = test.grid$y1, y2 = test.grid$y2,
  func = "gaussian", approx.method = "CC", boundary = boundary, max.prob = max.prob,
  initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs, params = params)

# Centroid-to-area transition matrix using Monte Carlo integration
CA.mc <- LatticeTransitionProbs(
  x1 = test.grid$x1, x2 = test.grid$x2, y1 = test.grid$y1, y2 = test.grid$y2,
  func = genexp.disp, approx.method = "CA", boundary = boundary, max.prob = max.prob,
  initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs, params = params)
# Centroid-to-area transition matrix using analytic results
CA.an <- LatticeTransitionProbs(
  x1 = test.grid$x1, x2 = test.grid$x2, y1 = test.grid$y1, y2 = test.grid$y2,
  func = "gaussian", approx.method = "CA", boundary = boundary, max.prob = max.prob,
  initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs, params = params)

# Area-to-centroid transition matrix using Monte Carlo integration
AC.mc <- LatticeTransitionProbs(
  x1 = test.grid$x1, x2 = test.grid$x2, y1 = test.grid$y1, y2 = test.grid$y2,
  func = genexp.disp, approx.method = "AC", boundary = boundary, max.prob = max.prob,
  initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs, params = params)
# Area-to-centroid transition matrix using analytic results
AC.an <- LatticeTransitionProbs(
  x1 = test.grid$x1, x2 = test.grid$x2, y1 = test.grid$y1, y2 = test.grid$y2,
  func = "gaussian", approx.method = "AC", boundary = boundary, max.prob = max.prob,
  initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs, params = params)

# Area-to-area transition matrix using Monte Carlo integration
AA.mc <- LatticeTransitionProbs(
  x1 = test.grid$x1, x2 = test.grid$x2, y1 = test.grid$y1, y2 = test.grid$y2,
  func = genexp.disp, approx.method = "AA", boundary = boundary, max.prob = max.prob,
  initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs, params = params)
# Area-to-area transition matrix using analytic results
AA.an <- LatticeTransitionProbs(
  x1 = test.grid$x1, x2 = test.grid$x2, y1 = test.grid$y1, y2 = test.grid$y2,
  func = "gaussian", approx.method = "AA", boundary = boundary, max.prob = max.prob,
  initial.step = initial.step, step.size = step.size, max.rel.var = max.rel.var, max.runs = max.runs, params = params)

# Calculate the differences between the analytic and Monte Carlo results
diff.CC <- CC.mc - CC.an
diff.CA <- CA.mc - CA.an
diff.AC <- AC.mc - AC.an
diff.AA <- AA.mc - AA.an

#=======================================
# https://www.r-bloggers.com/continuous-dispersal-on-a-discrete-lattice/

# General function to take in a lattice and disperse
## according to a user provided dispersal kernel
## Author: Corey Chivers
lat_disp<-function(pop,kernel,...)
{
  lattice_size<-dim(pop)
  new_pop<-array(0,dim=lattice_size)
  for(i in 1:lattice_size[1])
  {
    for(j in 1:lattice_size[2])
    {
      N<-pop[i,j]
      dist<-kernel(N,...)
      theta<-runif(N,0,2*pi)
      x<-cos(theta)*dist
      y<-sin(theta)*dist
      
      for(k in 1:N)
      {
        x_ind<-round(i+x[k]) %% lattice_size[1]
        y_ind<-round(j+y[k]) %% lattice_size[2]
        new_pop[x_ind,y_ind]<-new_pop[x_ind,y_ind]+1
      }
    }
  }
  return(new_pop)
}

############## Run and plot #######################

## Custom colour ramp
colours<-c('grey','blue','black')
cus_col<-colorRampPalette(colors=colours, space = c("rgb", "Lab"),interpolate = c("linear", "spline"))

## Initialize population array
Time=35
pop<-array(0,dim=c(Time,50,50))
pop[1,25,25]<-10000

### Normal Kernel ###
par(mfrow=c(1,1))
for(i in 2:Time)
{
  image(pop[i-1,,],col=cus_col(100),xaxt='n',yaxt='n')
  pop[i,,]<-lat_disp(pop[i-1,,],kernel=rnorm,mean=0,sd=1)
}

## Plot
png('normal_kern.png', width = 800, height = 800)
par(mfrow=c(2,2),pty='s',omi=c(0.1,0.1,0.5,0.1),mar=c(2,0,2,0))
times<-c(5,15,25,35)
for(i in times)
  image(pop[i-1,,],
        col=cus_col(100),
        xaxt='n',
        yaxt='n',
        useRaster=TRUE,
        main=paste("Time =",i))

mtext("Gaussian Kernel",outer=TRUE,cex=1.5)
dev.off()

### Exponential Kernel ###
par(mfrow=c(1,1))
for(i in 2:Time)
{
  image(pop[i-1,,],col=cus_col(100),xaxt='n',yaxt='n')
  pop[i,,]<-lat_disp(pop[i-1,,],kernel=rexp,rate=1)
}

## Plot
png('exp_kern.png',  width = 800, height = 800)
par(mfrow=c(2,2),pty='s',omi=c(0.1,0.1,0.5,0.1),mar=c(2,0,2,0))
times<-c(5,15,25,35)
for(i in times)
  image(pop[i-1,,],
        col=cus_col(100),
        xaxt='n',
        yaxt='n',
        useRaster=TRUE,
        main=paste("Time =",i))

mtext("Exponential Kernel",outer=TRUE,cex=1.5)
dev.off()




 