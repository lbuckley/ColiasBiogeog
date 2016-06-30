 #Evolutionary Model across biogeography

library(foreach)
 
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

#=====================================================
##  PLOT OUT
inds=1:137
cols=rainbow(841)

plot( years[inds],abs.mean[inds,1,1], type="l", ylim=range(0.52, 0.7) )

for(cell.k in 2:841){
 points( years[inds],abs.mean[inds,cell.k,1], type="l", col= cols[cell.k]) 
}

#------------------
#Lambda
plot( years[inds],lambda.mean[inds,1], type="l")

for(cell.k in 2:841){
  points( years[inds],lambda.mean[inds,cell.k], type="l", col= cols[cell.k]) 
}

DO:
#CALCULATE ABS SLOPE BY DECADE?
#ASSUME START AT OPTIMAL ABS. TOO COLD?
  
  