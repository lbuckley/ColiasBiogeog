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
#Load data
    
fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\"

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

#--------------------------
#Read optimal absorptivity
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")
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
#Load absorptivity and lambda values
#  abs.mean #dims: yr.k, cell.k, gen.k, scen.k:no plast, plast, only plast, metrics: abssample, absmid, rn, Babsmid, Brn)
#  lambda.mean= array(NA, dim=c(length(years),nrow(pts.sel), 3, 5)) #dims: yr.k, cell.k, gen.k, scen.k
# scen.k: plast0evol0, plast1evol0, plast0evol1, plast1evol1, plast1evol1rnevol1

abs.mean <- readRDS("absmean.abs")
lambda.mean <- readRDS("lambdamean.abs")

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
  
  # grid.newpage()
  grid.draw(combined)
  
  # return gtable invisibly
  invisible(combined)
  
}

#=====================================================
##  PLOTs

#Fig 1. elev vs year for Jadult, Tpup, Tad    
  
  #dim(pup.temps)= 12 150 841 3
inds=1:150
years= years[inds]

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
  
  #------------------
  
  #sample sites to faciliate visualization
  #s.inds=sort(base::sample(1:nrow(phen),50000))  
  phen1= phen #[s.inds,]
  
  #Interpolate
  s=interp(phen1$year,phen1$elev,phen1$Jadult, duplicate="strip")
  
  gdat <- interp2xyz(s, data.frame=TRUE)
  
  plot.Jad= ggplot(gdat) + 
    aes(x = x, y = y, z = z, fill = z) + 
    geom_tile() + 
    scale_fill_distiller(palette="Spectral", na.value="white", name="Jadult") +
    theme_bw(base_size=16)+xlab("year")+ylab("elevation (m)")+theme(legend.position="bottom")+ xlim(150, 270)
 
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
  #s.inds=sort(base::sample(1:nrow(phen),10000))  
  #phen1= phen[s.inds,]
  phen1= na.omit(phen)
  
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
  #s.inds=sort(base::sample(1:nrow(phen),50000))  
  #phen1= phen[s.inds,]
  
  phen1= na.omit(phen)
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
  
  #=====================================================
  
#Fig 2. PLOT FITNESS CURVES
  #Lambda[years, sites, abs, gen, metrics: Lambda, FAT,Egg Viability]
  
  abs= seq(0.4,0.7,0.05)
  
  for(gen.k in 1:2){
  
#group by elevation
  dat= Lambda[,,,gen.k,1]
  
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

  #-------------------
  #surface version of figure
  
  #drop NA lambdas
  fc= fc[which(!is.na(fc$lambda)),]
  
  fc$abs= abs[as.numeric(fc$abs)]  
  fc$year= as.numeric(fc$year)
  
  fc$yr.cut= cut(fc$elev, breaks=c(1950,1980,2010,2040,2070,2100),labels=c("1950_1980","2010","2040","2070","2100") )
  
  #pick time period
  #mean across years
  fc1= ddply(fc, .(ind,elev,yr.cut,abs), summarize, lambda=mean(lambda,na.rm=TRUE))
  fc2= fc1[which(fc1$yr.cut=="2010"),]
  
  #yrs= c(1950,1980,2010,2040,2070,2100)
  yrs= c(1980,2010,2040,2070)
  
  for(yr.k in 1:length(yrs)){
    
    #pick year
    fc2= fc[which(fc$year==yrs[yr.k] ),]
    
    #Interpolate
    s=interp(fc2$elev,fc2$abs,fc2$lambda, duplicate="mean")
    gdat <- interp2xyz(s, data.frame=TRUE)
    
    plot.lambda= ggplot(gdat) + 
      aes(x = x, y = y, z = z, fill = z) + 
      geom_tile() + 
      scale_fill_distiller(palette="Spectral", na.value="white", name="lambda") +
      theme_bw(base_size=16)+xlab("elevation (m)")+ylab("absorptivity")+theme(legend.position="bottom") +annotate("text", x=2000,y=0.68, label= paste(yrs[yr.k],"gen",gen.k,sep=" "), size=5)
    
    if(yr.k==1 & gen.k==1)plot.lambda1980.g1= plot.lambda
    if(yr.k==2 & gen.k==1)plot.lambda2010.g1= plot.lambda
    if(yr.k==3 & gen.k==1)plot.lambda2040.g1= plot.lambda
    if(yr.k==4 & gen.k==1)plot.lambda2070.g1= plot.lambda
    
    if(yr.k==1 & gen.k==2)plot.lambda1980.g2= plot.lambda
    if(yr.k==2 & gen.k==2)plot.lambda2010.g2= plot.lambda
    if(yr.k==3 & gen.k==2)plot.lambda2040.g2= plot.lambda
    if(yr.k==4 & gen.k==2)plot.lambda2070.g2= plot.lambda
    
  } #end loop years
  } #end loop generation
  
  #-------------------
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Fig2_FitnessCurves_elevs.pdf", height = 8, width = 12)
  
  grid_arrange_shared_legend(plot.lambda1980.g1,plot.lambda2010.g1,plot.lambda2040.g1,plot.lambda2070.g1,plot.lambda1980.g2,plot.lambda2010.g2,plot.lambda2040.g2,plot.lambda2070.g2, ncol = 4, nrow = 2)
  
  dev.off()  
  
#=======================================
#Fig 3. Fitness change by scenario
  
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
  l1$scen.name= scens[l1$scen+1]
  l2$scen.name= scens[l2$scen+1]
  l3$scen.name= scens[l3$scen+1]
  
  #-------------------
  #surface version of figure
  
  #LAMBDA
  for(gen.k in 1:2){
    
    if(gen.k==1) fcl=l1
    if(gen.k==2) fcl=l2
    
    for(scen.k in 1:4){
      
      #pick scenarios
      fc= fcl[which(fcl$scen==scen.k),]
      
      #drop outliers
      fc= fc[which(fc$lambda> (-1) & fc$lambda< 1),]
      
      #drop NA lambdas
      fc= fc[which(!is.na(fc$lambda)),]
      
      #sample
      if(nrow(fc)>50000 ) inds= base::sample(1:nrow(fc),50000, replace=FALSE)
      fc=fc[inds,]
      
      #Interpolate
      s=interp(fc$year,fc$elev,fc$lambda, duplicate="mean")
      gdat <- interp2xyz(s, data.frame=TRUE)
      
      plot.lambda= ggplot(gdat) + 
        aes(x = x, y = y, z = z, fill = z) + 
        geom_tile() + 
        scale_fill_distiller(palette="Spectral", na.value="white", name="lambda") +
        theme_bw(base_size=16)+xlab("year")+ylab("elevation (m)")+theme(legend.position="bottom") +annotate("text", x=2000,y=1500, label= fc$scen.name[1], size=5)
      
      if(scen.k==1 & gen.k==1)plot.lambda.scen1.g1= plot.lambda
      if(scen.k==2 & gen.k==1)plot.lambda.scen2.g1= plot.lambda
      if(scen.k==3 & gen.k==1)plot.lambda.scen3.g1= plot.lambda
      if(scen.k==4 & gen.k==1)plot.lambda.scen4.g1= plot.lambda
      
      if(scen.k==1 & gen.k==2)plot.lambda.scen1.g2= plot.lambda
      if(scen.k==2 & gen.k==2)plot.lambda.scen2.g2= plot.lambda
      if(scen.k==3 & gen.k==2)plot.lambda.scen3.g2= plot.lambda
      if(scen.k==4 & gen.k==2)plot.lambda.scen4.g2= plot.lambda
      
    } #end loop scenario
    
  } #end loop generation
  
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
  l1$scen.name= scens[l1$scen+1]
  l2$scen.name= scens[l2$scen+1]
  l3$scen.name= scens[l3$scen+1]
  
  #-------------------
  #surface version of figure
  
  #LAMBDA
  for(gen.k in 1:1){
    
    if(gen.k==1) fcl=l1
    if(gen.k==2) fcl=l2
    
    for(scen.k in 1:4){
      
      #pick scenarios
      fc= fcl[which(fcl$scen==scen.k),]
      
      #drop outliers
      fc= fc[which(fc$abs> (-0.2) & fc$abs< 0.2),]
      
      #drop NA lambdas
      fc= fc[which(!is.na(fc$abs)),]
      
      #sample
      if(nrow(fc)>50000 ) inds= base::sample(1:nrow(fc),50000, replace=FALSE)
      fc=fc[inds,]
      
      #Interpolate
      s=interp(fc$year,fc$elev,fc$abs, duplicate="mean")
      gdat <- interp2xyz(s, data.frame=TRUE)
      
      plot.abs= ggplot(gdat) + 
        aes(x = x, y = y, z = z, fill = z) + 
        geom_tile() + 
        scale_fill_distiller(palette="Spectral", na.value="white", name="lambda") +
        theme_bw(base_size=16)+xlab("year")+ylab("elevation (m)")+theme(legend.position="bottom") +annotate("text", x=2000,y=1500, label= fc$scen.name[1], size=5)
      
      if(scen.k==1 & gen.k==1)plot.abs.scen1.g1= plot.abs
      if(scen.k==2 & gen.k==1)plot.abs.scen2.g1= plot.abs
      if(scen.k==3 & gen.k==1)plot.abs.scen3.g1= plot.abs
      if(scen.k==4 & gen.k==1)plot.abs.scen4.g1= plot.abs
      
      if(scen.k==1 & gen.k==2)plot.abs.scen1.g2= plot.abs
      if(scen.k==2 & gen.k==2)plot.abs.scen2.g2= plot.abs
      if(scen.k==3 & gen.k==2)plot.abs.scen3.g2= plot.abs
      if(scen.k==4 & gen.k==2)plot.abs.scen4.g2= plot.abs
      
    } #end loop scenario
    
  } #end loop generation
  
  #-------------------
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Fig3_LambdaSurf_scen.pdf", height = 8, width = 12)
  
  grid_arrange_shared_legend(plot.lambda.scen1.g1,plot.lambda.scen2.g1,plot.lambda.scen3.g1,plot.lambda.scen4.g1,plot.lambda.scen1.g2,plot.lambda.scen2.g2,plot.lambda.scen3.g2,plot.lambda.scen4.g2, ncol = 4, nrow = 2)
  
  dev.off()  
  
  #-----------------
  setwd(paste(fdir,"figures\\",sep="") )
  pdf("Fig3_AbsSurf_scen.pdf", height = 4, width = 9)
  
  grid_arrange_shared_legend(plot.abs.scen2.g1, plot.abs.scen3.g1, plot.abs.scen4.g1, ncol = 3, nrow = 1)
  
  dev.off() 
  
#===============================================
#MAPS LAMBDA AND ABS RELATIVE TO BASELINE (2010?)

#CURRENTLY na.rm=TRUE, change?

#LAMBDA
#average lambdas across time period
#first generation
lambda.diff= lambda.mean[,,1,]

#change erroneous values
lambda.diff[which(lambda.diff< (-10))]=NA
lambda.diff[which(lambda.diff> 10)]=NA

per1= colMeans(lambda.diff[which(years %in% 1951:1980),,], na.rm = TRUE, dims = 1)
per2= colMeans(lambda.diff[which(years %in% 1981:2010),,], na.rm = TRUE, dims = 1)
per3= colMeans(lambda.diff[which(years %in% 2011:2040),,], na.rm = TRUE, dims = 1)
per4= colMeans(lambda.diff[which(years %in% 2041:2070),,], na.rm = TRUE, dims = 1)

#translate to difference from no plasticity or evolution
lper1s= per1 -per1[,1]
lper2s= per2 -per2[,1]
lper3s= per3 -per3[,1]
lper4s= per4 -per4[,1]

##fixed tiem period
#lper1s= per1 -per1[,1]
#lper2s= per2 -per1[,1]
#lper3s= per3 -per1[,1]
#lper4s= per4 -per1[,1]

#----------
#second generation

lambda.diff= lambda.mean[,,2,]
per1= colMeans(lambda.diff[which(years %in% 1951:1980),,], na.rm = TRUE, dims = 1)
per2= colMeans(lambda.diff[which(years %in% 1981:2010),,], na.rm = TRUE, dims = 1)
per3= colMeans(lambda.diff[which(years %in% 2011:2040),,], na.rm = TRUE, dims = 1)
per4= colMeans(lambda.diff[which(years %in% 2041:2070),,], na.rm = TRUE, dims = 1)

#translate to difference from no plasticity or evolution
lper1s.gen2= per1 -per1[,1]
lper2s.gen2= per2 -per2[,1]
lper3s.gen2= per3 -per3[,1]
lper4s.gen2= per4 -per4[,1]

##fixed time period
#lper1s.gen2= per1 -per1[,1]
#lper2s.gen2= per2 -per1[,1]
#lper3s.gen2= per3 -per1[,1]
#lper4s.gen2= per4 -per1[,1]

#correct NaNs
lper1s[is.nan(lper1s)] <- NA
lper2s[is.nan(lper2s)] <- NA
lper3s[is.nan(lper3s)] <- NA
lper4s[is.nan(lper4s)] <- NA
#----------

#ABS
#average abs across time period
#first generation
lambda.diff= abs.mean[,,1,,2]#abs mid

#change erroneous values
lambda.diff[which(lambda.diff>0.8 | lambda.diff<0.2)]=NA

per1= colMeans(lambda.diff[which(years %in% 1951:1980),,], na.rm = TRUE, dims = 1)
per2= colMeans(lambda.diff[which(years %in% 1981:2010),,], na.rm = TRUE, dims = 1)
per3= colMeans(lambda.diff[which(years %in% 2011:2040),,], na.rm = TRUE, dims = 1)
per4= colMeans(lambda.diff[which(years %in% 2041:2070),,], na.rm = TRUE, dims = 1)

#translate to difference from 1981-2010 for no plasticity or evolution
aper1s= per1 -per1[,1]
aper2s= per2 -per2[,1]
aper3s= per3 -per3[,1]
aper4s= per4 -per4[,1]
#correct NaNs
aper1s[is.nan(aper1s)] <- NA
aper2s[is.nan(aper2s)] <- NA
aper3s[is.nan(aper3s)] <- NA
aper4s[is.nan(aper4s)] <- NA

#determine breaks
l.break= rbind(lper1s, lper2s, lper3s, lper4s)
l.breaks=quantile(l.break, probs=seq(0,1,0.1), na.rm=TRUE)

a.break= rbind(aper1s, aper2s, aper3s, aper4s)
a.breaks=quantile(a.break, probs=seq(0,1,0.1), na.rm=TRUE)

l.breaks= unique(l.breaks)
a.breaks= unique(a.breaks)

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
  
  #bin
  lper1$lambda.bins= cut(lper1$lambda, breaks=l.breaks)
  lper2$lambda.bins= cut(lper2$lambda, breaks=l.breaks)
  lper3$lambda.bins= cut(lper3$lambda, breaks=l.breaks)
  lper4$lambda.bins= cut(lper4$lambda, breaks=l.breaks)
  
  #set up map
  bbox <- ggmap::make_bbox(lon, lat, lper1, f = 0.1)
  map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')
  map1=ggmap(map_loc, margins=FALSE)
  
  #OLD VERSION
  #lper1.map<- map1 + geom_raster(data=lper1, aes(fill = lambda), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=l.breaks, labels=l.lab )+ coord_fixed() + theme(legend.position="right")

  lper1.map<- map1 + geom_raster(data=lper1, aes(fill = lambda.bins), alpha=0.5)+ coord_cartesian()+ scale_fill_brewer(name="lambda", palette="PuOr")+ coord_fixed() + theme(legend.position="right")

  lper2.map<- map1 + geom_raster(data=lper2, aes(fill = lambda.bins), alpha=0.5)+ coord_cartesian()+ scale_fill_brewer(name="lambda", palette="PuOr")+ coord_fixed() + theme(legend.position="right")
  
  lper3.map<- map1 + geom_raster(data=lper3, aes(fill = lambda.bins), alpha=0.5)+ coord_cartesian()+ scale_fill_brewer(name="lambda", palette="PuOr")+ coord_fixed() + theme(legend.position="right")
  
  lper4.map<- map1 + geom_raster(data=lper4, aes(fill = lambda.bins), alpha=0.5)+ coord_cartesian()+ scale_fill_brewer(name="lambda", palette="PuOr")+ coord_fixed() + theme(legend.position="right")
  
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
  
  #bin
  aper1$abs.bins= cut(aper1$abs, breaks=a.breaks)
  aper2$abs.bins= cut(aper2$abs, breaks=a.breaks)
  aper3$abs.bins= cut(aper3$abs, breaks=a.breaks)
  aper4$abs.bins= cut(aper4$abs, breaks=a.breaks)
  
  #set up map
  bbox <- ggmap::make_bbox(lon, lat, aper1, f = 0.1)
  map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')
  map1=ggmap(map_loc, margins=FALSE)

#OLD VERSION  
#  aper1.map<- map1 + geom_raster(data=aper1, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=a.breaks, labels=a.lab)+ coord_fixed() + theme(legend.position="right")
  
  aper1.map<- map1 + geom_raster(data=aper1, aes(fill = abs.bins), alpha=0.5)+ coord_cartesian()+ scale_fill_brewer(name="abs", palette="PuOr")+ coord_fixed() + theme(legend.position="right")
  
  aper2.map<- map1 + geom_raster(data=aper2, aes(fill = abs.bins), alpha=0.5)+ coord_cartesian()+ scale_fill_brewer(name="abs", palette="PuOr")+ coord_fixed() + theme(legend.position="right")
  
  aper3.map<- map1 + geom_raster(data=aper3, aes(fill = abs.bins), alpha=0.5)+ coord_cartesian()+ scale_fill_brewer(name="abs", palette="PuOr")+ coord_fixed() + theme(legend.position="right")
  
  aper4.map<- map1 + geom_raster(data=aper4, aes(fill = abs.bins), alpha=0.5)+ coord_cartesian()+ scale_fill_brewer(name="abs", palette="PuOr")+ coord_fixed() + theme(legend.position="right")
  
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
#quilt.plot(lper1$lon,lper1$lat, lper1$lambda)

#Try interpolation
aper2.map<- map1 + geom_raster(data=aper2, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ scale_fill_gradientn(colours=matlab.like(10), breaks=a.breaks, labels=a.lab)+ coord_fixed() + theme(legend.position="right")

aper2.map<- map1 + geom_raster(data=aper2, interpolate=TRUE, aes(fill = abs), alpha=0.5)+ coord_cartesian()

#----------------------------

#FIG. 4: ABS MAP
setwd(paste(fdir,"figures\\",sep="") )
pdf("Fig4_Abs_map.pdf", height = 4, width = 8)

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
#FIG. 5: LAMBDA MAP

pdf("Fig5_Lambda_map.pdf", height = 4, width = 8)

grid_arrange_shared_legend(lper1.3,lper2.3,lper3.3,lper4.3,lper1.4,lper2.4,lper3.4,lper4.4,lper1.5,lper2.5,lper3.5,lper4.5, ncol = 4, nrow = 3)

dev.off()

#lper1.1,lper2.1,lper3.1,lper4.1,lper1.2,lper2.2,lper3.2,lper4.2, 


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
