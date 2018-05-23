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
library(RColorBrewer)

#pick projection
proj.k=2
projs=c("bcc-csm","ccsm4","gfdl")

years= 1950:2099
#-----------------------------
#Load data
    
#fdir= "C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\"
fdir= "/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/" 

#Read points
#setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")
setwd("/Volumes/GoogleDrive/My\ Drive/Buckley/Work/ColiasBiogeog/OUT/")
pts.sel= read.csv( paste("COpoints.csv", sep="") ) #_",projs[proj.k],"
  
#Read lambdas and pupal temps
#setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")
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
yr.inds=1:150

#--------------------------
#Read optimal absorptivity
#setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")
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

inds=1:150
years= years[inds]

#Fig 0. Maps
#a: color= elevation
#b: color= Tpupal (1st gen)
#c: color= DOY adult (1st gen)

#Jadult
phen1= cbind(pts.sel, t(pup.temps["Jadult",yr.inds,,1]) ) 
phen1$gen=1
#extract 1950-1980
phen1=phen1[,1:(9+31)]
#mean across years
phen1$Jadult= rowMeans(phen1[,9:40] )

#Tpupal
phen2= cbind(pts.sel, t(pup.temps["Tpup",yr.inds,,1]) ) 
phen2$gen=1
#extract 1950-1980
phen2=phen2[,1:(9+31)]
#mean across years
phen2$Tpupal= rowMeans(phen2[,9:40] )

#Jadult and Tpupal
phen= cbind(phen1[,c(1:8)], phen1$Jadult, phen2$Tpupal)

#-----------------

#set up map
bbox <- ggmap::make_bbox(lon, lat, phen, f = 0.1)
map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')
map1=ggmap(map_loc, margins=FALSE)

#elevation
Elev<- map1 + geom_raster(data=phen, aes(fill = elev), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="bottom")+xlab("")+ylab("Latitude (°)")+ scale_fill_gradientn(colours=matlab.like(10), name="Elevation (m)")+ theme(legend.key.width=unit(0.8,"cm"))

#Jadult
Jadult<- map1 + geom_raster(data=phen, aes(fill = phen1$Jadult), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="bottom")+xlab("Longitude (°)")+ylab(NULL)+ scale_fill_gradientn(colours=matlab.like(10), name="Phenology (doy)")+ theme(legend.key.width=unit(0.8,"cm"))

#Tpupal
Tpupal<- map1 + geom_raster(data=phen, aes(fill = phen2$Tpupal), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="bottom")+xlab("")+ylab(NULL)+ scale_fill_gradientn(colours=matlab.like(10), name="Pupal temperature (°C)")+ theme(legend.key.width=unit(0.8,"cm"))

#plot together
library(cowplot)

setwd(paste(fdir, "figures/", sep=""))
pdf("Fig0_map.pdf", height=8, width=12)

plot_grid(Elev, Jadult, Tpupal, align = "h", ncol = 3, rel_widths = c(1.1,1,1))

dev.off()

#---------------------------------

#Fig 1. elev vs year for Jadult, Tpup, Tad    

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

#JADULT
#Jadult by gen
phen1= cbind(pts.sel, t(pup.temps["Jadult",yr.inds,,1]) ) 
phen1$gen=1
phen2= cbind(pts.sel, t(pup.temps["Jadult",yr.inds,,2]) ) 
phen2$gen=2
phen3= cbind(pts.sel, t(pup.temps["Jadult",yr.inds,,3]) ) 
phen3$gen=3
phen= rbind(phen1, phen2, phen3)

phen= gather(phen, "year", "Jadult",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]

#replace Jadult=273 with NA
phen$Jadult[which(phen$Jadult==273)]=NA
#remove NAs
#phen= phen[which(!is.na(phen$Jadult)),]

#time periods
phen$yr.cut= cut(phen$year, breaks=c(1950,1980,2010,2040,2070,2100),labels=c("1950_1980","2010","2040","2070","2100") )

#restrict to sites with 3 generations in all years
phen.3gen= ddply(phen, .(ind,elev), summarize, Jadult=mean(Jadult,na.rm=FALSE))
phen.3gen= na.omit(phen.3gen)
phen1= phen[which(phen$ind %in% phen.3gen$ind),]

#mean across years
phen2= ddply(phen1, .(ind,elev, yr.cut,gen), summarize, Jadult=mean(Jadult,na.rm=TRUE))
#remove NAs
phen2= phen2[which(!is.na(phen2$Jadult)),]

#select time
phen1980= phen2[which(phen2$yr.cut=="1950_1980"),]
phen2040= phen2[which(phen2$yr.cut=="2040"),]

##Save Tpupal
Jad1980= phen1980

#match
phen1980$id= paste(phen1980$ind, phen1980$gen)
phen2040$id= paste(phen2040$ind, phen2040$gen)
match1= match(phen1980$id, phen2040$id)
phen2040$Jadult1980= phen1980$Jadult[match1]
#diff
phen2040$Jadult_diff= phen2040$Jadult - phen2040$Jadult1980

#----------------------------
#PLOT

#Interpolate
s=interp(x=phen1980$gen,y=phen1980$elev,z=phen1980$Jadult, duplicate="mean", xo=c(1,2,3))

gdat <- interp2xyz(s, data.frame=TRUE)

plot.Jad= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="doy") +
  theme_classic(base_size=16)+xlab(NULL)+ylab("elevation (m)")+theme(legend.position="bottom")+labs(title= "baseline (1951-1980)")+
  theme(plot.title = element_text(size = 16))+
  theme(legend.margin=margin(-6,0,0,0) )+ theme(legend.key.width=unit(0.8,"cm")) +theme(axis.title.x=element_blank())

#+annotate("text", x=1,y=3000, label= "1951-1980", size=5)

#----------------
s=interp(x=phen2040$gen,y=phen2040$elev,z=phen2040$Jadult_diff, duplicate="mean", xo=c(1,2,3))

gdat <- interp2xyz(s, data.frame=TRUE)

plot.Jad2040= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="doy") +
  theme_classic(base_size=16)+xlab("")+ylab("elevation (m)")+theme(legend.position="bottom")+labs(title= expression(paste(Delta, " from baseline")) )+
  theme(plot.title = element_text(size = 16))

#-------------------------------
#Tpup
#Tpup by gen
phen1= cbind(pts.sel, t(pup.temps["Tpup",yr.inds,,1]) ) 
phen1$gen=1
phen2= cbind(pts.sel, t(pup.temps["Tpup",yr.inds,,2]) ) 
phen2$gen=2
phen3= cbind(pts.sel, t(pup.temps["Tpup",yr.inds,,3]) ) 
phen3$gen=3
phen= rbind(phen1, phen2, phen3)

phen= gather(phen, "year", "Tpup",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]

#remove NAs
#phen= phen[which(!is.na(phen$Tpup)),]

#time periods
phen$yr.cut= cut(phen$year, breaks=c(1950,1980,2010,2040,2070,2100),labels=c("1950_1980","2010","2040","2070","2100") )

#restrict to sites with 3 generations in all years
phen.3gen= ddply(phen, .(ind,elev), summarize, Tpup=mean(Tpup,na.rm=FALSE))
phen.3gen= na.omit(phen.3gen)
phen1= phen[which(phen$ind %in% phen.3gen$ind),]

#mean across years
phen2= ddply(phen1, .(ind,elev, yr.cut,gen), summarize, Tpup=mean(Tpup,na.rm=TRUE))

#select time
phen1980= phen2[which(phen2$yr.cut=="1950_1980"),]
phen2040= phen2[which(phen2$yr.cut=="2040"),]

##Save Tpupal
Tpup1980= phen1980

#match
phen1980$id= paste(phen1980$ind, phen1980$gen)
phen2040$id= paste(phen2040$ind, phen2040$gen)
match1= match(phen1980$id, phen2040$id)
phen2040$Tpup1980= phen1980$Tpup[match1]
#diff
phen2040$Tpup_diff= phen2040$Tpup - phen2040$Tpup1980

#Interpolate
s=interp(x=phen1980$gen,y=phen1980$elev,z=phen1980$Tpup, duplicate="mean", xo=c(1,2,3))

gdat <- interp2xyz(s, data.frame=TRUE)

plot.Tpup= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="Tpupal (°C)") +
  theme_classic(base_size=16)+xlab(NULL)+ylab(NULL)+theme(legend.position="bottom")+labs(title= "")+
  theme(legend.margin=margin(-6,0,0,0) )+theme(axis.title.x=element_blank())

#----------------
#remove anomolous sites with large temp declines to prevent legend issues
s=interp(x=phen2040$gen,y=phen2040$elev,z=phen2040$Tpup_diff, duplicate="mean", xo=c(1,2,3))

gdat <- interp2xyz(s, data.frame=TRUE)

plot.Tpup2040= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="Tpupal (°C)", limits=c(-1.11,2.6) ) +
  theme_classic(base_size=16)+xlab("generation")+ylab(NULL)+theme(legend.position="bottom")+labs(title= "")
#set limit to omit 1 outlier

#-------------------------------
#Tad
#Tad by gen
phen1= cbind(pts.sel, t(pup.temps["Tad",yr.inds,,1]) ) 
phen1$gen=1
phen2= cbind(pts.sel, t(pup.temps["Tad",yr.inds,,2]) ) 
phen2$gen=2
phen3= cbind(pts.sel, t(pup.temps["Tad",yr.inds,,3]) ) 
phen3$gen=3
phen= rbind(phen1, phen2, phen3)

phen= gather(phen, "year", "Tad",9:(length(years)+8))
phen$year= years[as.numeric(phen$year)]

#remove NAs
#phen= phen[which(!is.na(phen$Tad)),]

#time periods
phen$yr.cut= cut(phen$year, breaks=c(1950,1980,2010,2040,2070,2100),labels=c("1950_1980","2010","2040","2070","2100") )

#restrict to sites with 3 generations in all years
phen.3gen= ddply(phen, .(ind,elev), summarize, Tad=mean(Tad,na.rm=FALSE))
phen.3gen= na.omit(phen.3gen)
phen1= phen[which(phen$ind %in% phen.3gen$ind),]

#mean across years
phen2= ddply(phen1, .(ind,elev, yr.cut,gen), summarize, Tad=mean(Tad,na.rm=TRUE))

#select time
phen1980= phen2[which(phen2$yr.cut=="1950_1980"),]
phen2040= phen2[which(phen2$yr.cut=="2040"),]

#match
phen1980$id= paste(phen1980$ind, phen1980$gen)
phen2040$id= paste(phen2040$ind, phen2040$gen)
match1= match(phen1980$id, phen2040$id)
phen2040$Tad1980= phen1980$Tad[match1]
#diff
phen2040$Tad_diff= phen2040$Tad - phen2040$Tad1980

#Interpolate
s=interp(x=phen1980$gen,y=phen1980$elev,z=phen1980$Tad, duplicate="mean", xo=c(1,2,3))

gdat <- interp2xyz(s, data.frame=TRUE)

plot.Tad= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="Tadult (°C)") +
  theme_classic(base_size=16)+xlab(NULL)+ylab(NULL)+theme(legend.position="bottom")+labs(title= "")+
  theme(legend.margin=margin(-6,0,0,0) )+theme(axis.title.x=element_blank())
#----------------
s=interp(x=phen2040$gen,y=phen2040$elev,z=phen2040$Tad_diff, duplicate="mean", xo=c(1,2,3))

gdat <- interp2xyz(s, data.frame=TRUE)

plot.Tad2040= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="Tadult (°C)") +
  theme_classic(base_size=16)+xlab("")+ylab(NULL)+theme(legend.position="bottom")+labs(title= "")

#=============================

# plot together
setwd(paste(fdir,"figures/", sep=""))

pdf("Fig1_FigJadTpupTad.pdf", height=8, width=12)

plot_grid(plot.Jad, plot.Tpup, plot.Tad, plot.Jad2040, plot.Tpup2040, plot.Tad2040, align = "h", nrow = 2, rel_heights = c(1,1), rel_widths = c(1.2,1,1))

dev.off()

  #=====================================================
  
#Fig 2. PLOT FITNESS CURVES
  #Lambda[years, sites, abs, gen, metrics: Lambda, FAT,Egg Viability]
  
  abs= seq(0.4,0.7,0.05)
  
#get legend
g_legend <-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend) }

#legend limits
dat= Lambda[,,,,1]
lambda.lims= range(dat, na.rm=TRUE)

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
  
  fc$yr.cut= cut(fc$year, breaks=c(1950,1980,2010,2040,2070,2100),labels=c("1950_1980","2010","2040","2070","2100") )
  
  #pick time period
  #mean across years
  fc1= ddply(fc, .(ind,elev,yr.cut,abs), summarize, lambda=mean(lambda,na.rm=TRUE))
  
  #yrs= c(1950,1980,2010,2040,2070,2100)
  yrs= c(1980,2010,2040,2070)
  yr.cuts= c("1950_1980","2010","2040","2070","2100")
  yr.labs= c("1950-1980","1981-2010","2011-2040","2041-2070","2071-2100")
  
  for(yr.k in 1:length(yrs)){
    
    #pick year
    fc2=  fc1[which(fc1$yr.cut==yr.cuts[yr.k] ),]
    #fc2= fc[which(fc$year==yrs[yr.k] ),]
    
    #Interpolate
    s=interp(fc2$elev,fc2$abs,fc2$lambda, duplicate="mean")
    gdat <- interp2xyz(s, data.frame=TRUE)
    
    #get legend
    if(yr.k==1 & gen.k==1){
      plot.lambda= ggplot(gdat) + 
        aes(x = x, y = y, z = z, fill = z) + 
        geom_tile() + 
        scale_fill_distiller(palette="Spectral", na.value="white", name="lambda", limits=lambda.lims) +theme(legend.position="right")
      
      leg = g_legend(plot.lambda)
}
          
    plot.lambda= ggplot(gdat) + 
      aes(x = x, y = y, z = z, fill = z) + 
      geom_tile() + 
      scale_fill_distiller(palette="Spectral", na.value="white", name="lambda", limits=lambda.lims) +
      theme_classic(base_size=16)+xlab(NULL)+ylab(NULL)+theme(legend.position="none") +labs(title= paste("generation", gen.k, "\n",yr.labs[yr.k],sep=" "), size=10)
    
    if(yr.k>1) plot.lambda= plot.lambda  +labs(title= paste("","\n",yr.labs[yr.k],sep=" "), size=10)
    
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
  
# +xlab("elevation (m)")+ylab("absorptivity")

  #-------------------
#lay out plots

f1 = plot_grid(plot.lambda1980.g1,plot.lambda2010.g1,plot.lambda2040.g1,plot.lambda2070.g1,plot.lambda1980.g2,plot.lambda2010.g2,plot.lambda2040.g2,plot.lambda2070.g2, ncol = 4, align="v")

f2= grid.arrange(
  arrangeGrob(f1,  
              bottom=grid::textGrob(label= "elevation (m)", gp= gpar(fontsize=24, col="black")),
              left=grid::textGrob(label= "absorptivity", rot=90, gp= gpar(fontsize=24, col="black"))),
  leg, 
  widths=c(9,1))

#---
  setwd(paste(fdir,"figures/",sep="") )
  pdf("Fig2_FitnessCurves_elevs.pdf", height = 7, width = 12)
  grid.draw(grobTree(rectGrob(gp=gpar(fill="white", lwd=0)), 
                     f2))
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
  
  #legend limits
  lambda.lims= c(-1,1)
  #lambda breaks
  lambda.breaks=quantile(lambda.diff, probs = seq(0, 1, length.out=7), na.rm=TRUE )
  
  scen.labs=  c("plast only", "evol only", "plast + evol", "+ evol of plast")
  
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
      
      #get legend
      if(scen.k==1 & gen.k==1){
        plot.lambda= ggplot(gdat) + 
          aes(x = x, y = y, z = z, fill = z) + 
          geom_tile() + theme(legend.position="right")+scale_fill_gradientn(colours=brewer.pal(7, "Spectral"),values=rescale(lambda.breaks), guide="colorbar", na.value="white",name="lambda", limits=lambda.lims)
        
        leg = g_legend(plot.lambda)
      }
      
      plot.lambda= ggplot(gdat) + 
        aes(x = x, y = y, z = z, fill = z) + 
        geom_tile() +scale_fill_gradientn(colours=brewer.pal(7, "Spectral"),values=rescale(lambda.breaks), guide="colorbar", na.value="white",name="lambda", limits=lambda.lims) +
        theme_classic(base_size=16)+xlab(NULL)+ylab(NULL)+theme(legend.position="none") +labs(title= paste(scen.labs[scen.k], "\n","generation", gen.k,sep=" "), size=10)
      
      if(gen.k>1) plot.lambda= plot.lambda  +labs(title= paste("generation", gen.k,sep=" "), size=10)
  
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
  
  #-------------------------
  #lay out plots
  
  f1 = plot_grid(plot.lambda.scen1.g1,plot.lambda.scen2.g1,plot.lambda.scen3.g1,plot.lambda.scen4.g1,plot.lambda.scen1.g2,plot.lambda.scen2.g2,plot.lambda.scen3.g2,plot.lambda.scen4.g2, ncol = 4, align="v")
  
  lambda.plot= grid.arrange(
    arrangeGrob(f1,  
                bottom=grid::textGrob(label= "year", gp= gpar(fontsize=24, col="black")),
                left=grid::textGrob(label= "elevation (m)", rot=90, gp= gpar(fontsize=24, col="black"))),
    leg, 
    widths=c(9,1))
  
  #====================================
#ABSORPTIVITIES
  
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
  
  #abs breaks
  lambda.diff2<- lambda.diff
  lambda.diff2[which(abs(lambda.diff2)>0.2) ]<-NA
  abs.breaks=quantile(lambda.diff2, probs = seq(0, 1, length.out=7), na.rm=TRUE )
  #abs limits
  abs.lims= range(lambda.diff2, na.rm=TRUE)
  
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
  
  scen.labs=  c("plast only", "evol only", "plast + evol", "+ evol of plast")
  
  #-------------------
  #surface version of figure
  
  #ABS
  for(gen.k in 1:2){
    
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
     
      #get legend
      if(scen.k==1 & gen.k==1){
        plot.abs= ggplot(gdat) + 
          aes(x = x, y = y, z = z, fill = z) + 
          geom_tile() + 
          scale_fill_gradientn(colours=brewer.pal(7, "Spectral"),values=rescale(abs.breaks), guide="colorbar", na.value="white",name="absorptivity", limits=abs.lims) +theme(legend.position="right")
        
        leg = g_legend(plot.abs)
      }
      
      plot.abs= ggplot(gdat) + 
        aes(x = x, y = y, z = z, fill = z) + 
        geom_tile() + 
        scale_fill_gradientn(colours=brewer.pal(7, "Spectral"),values=rescale(abs.breaks), guide="colorbar", na.value="white",name="absorptivity", limits=abs.lims) +
        theme_classic(base_size=16)+xlab(NULL)+ylab(NULL)+theme(legend.position="none") +labs(title= scen.labs[scen.k], size=10)

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
  
  #-----------------
  #Absorptivity
  #lay out plots
  
 ## check generation effect 
 # plot_grid(plot.abs.scen4.g1, plot.abs.scen4.g2)
  
  f1 = plot_grid(frame(),plot.abs.scen2.g1, plot.abs.scen3.g1, plot.abs.scen4.g1, ncol = 4, align="v")
  
  abs.plot= grid.arrange(
    arrangeGrob(f1,  
                bottom=grid::textGrob(label= "year", gp= gpar(fontsize=24, col="black")),
                left=grid::textGrob(label= "elevation (m)", rot=90, gp= gpar(fontsize=24, col="black"))),
    leg, 
    widths=c(9,1))
  
   #-------------------
 
  setwd(paste(fdir,"figures/",sep="") )
  pdf("Fig3.pdf", height = 8, width = 12)
  grid.draw(grobTree(rectGrob(gp=gpar(fill="white", lwd=0)), 
                     grid.arrange(lambda.plot, abs.plot, ncol=1, heights=c(2,1.2) )))
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

scen.labs= c("none", "plast only","evol only","plast + evol","+ evol of plast")

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
  
  #first map
  lper1.map<- map1 + geom_raster(data=lper1, aes(fill = lambda), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="right")+xlab(NULL)+ylab(NULL)+labs(title= scen.labs[scen.k], size=10)+
    scale_fill_gradientn(colours=brewer.pal(8, "RdYlBu"),values=rescale(l.breaks), guide="colorbar", na.value="white",name="lambda", limits=range(l.breaks))
  
  #get legend
  l.leg= get_legend(lper1.map)
  
  lper1.map<- lper1.map + theme(legend.position="none")
  
  #next maps
  lper2.map<- map1 + geom_raster(data=lper2, aes(fill = lambda), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="none")+xlab(NULL)+ylab(NULL)+labs(title= scen.labs[scen.k], size=10)+
    scale_fill_gradientn(colours=brewer.pal(8, "RdYlBu"),values=rescale(l.breaks), guide="colorbar", na.value="white",name="lambda", limits=range(l.breaks))
  
  lper3.map<- map1 + geom_raster(data=lper3, aes(fill = lambda), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="none")+xlab(NULL)+ylab(NULL)+labs(title= scen.labs[scen.k], size=10)+
    scale_fill_gradientn(colours=brewer.pal(8, "RdYlBu"),values=rescale(l.breaks), guide="colorbar", na.value="white",name="lambda", limits=range(l.breaks))
  
  lper4.map<- map1 + geom_raster(data=lper4, aes(fill = lambda), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="none")+xlab(NULL)+ylab(NULL)+labs(title= scen.labs[scen.k], size=10)+
    scale_fill_gradientn(colours=brewer.pal(8, "RdYlBu"),values=rescale(l.breaks), guide="colorbar", na.value="white",name="lambda", limits=range(l.breaks))

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

  #first map
  aper1.map<- map1 + geom_raster(data=aper1, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="right")+xlab(NULL)+ylab(NULL)+labs(title= scen.labs[scen.k], size=10)+
    scale_fill_gradientn(colours=brewer.pal(8, "RdYlBu"),values=rescale(a.breaks), guide="colorbar", na.value="white",name="absorptivity", limits=range(a.breaks))
  
  #get legend
  a.leg= get_legend(aper1.map)
  
  aper1.map<- aper1.map + theme(legend.position="none")
  
  #next maps
  aper2.map<- map1 + geom_raster(data=aper2, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="none")+xlab(NULL)+ylab(NULL)+labs(title= scen.labs[scen.k], size=10)+
    scale_fill_gradientn(colours=brewer.pal(8, "RdYlBu"),values=rescale(a.breaks), guide="colorbar", na.value="white",name="absorptivity", limits=range(a.breaks))
  
  aper3.map<- map1 + geom_raster(data=aper3, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="none")+xlab(NULL)+ylab(NULL)+labs(title= scen.labs[scen.k], size=10)+
    scale_fill_gradientn(colours=brewer.pal(8, "RdYlBu"),values=rescale(a.breaks), guide="colorbar", na.value="white",name="absorptivity", limits=range(a.breaks))
  
  aper4.map<- map1 + geom_raster(data=aper4, aes(fill = abs), alpha=0.5)+ coord_cartesian()+ coord_fixed() + theme(legend.position="none")+xlab(NULL)+ylab(NULL)+labs(title= scen.labs[scen.k], size=10)+
    scale_fill_gradientn(colours=brewer.pal(8, "RdYlBu"),values=rescale(a.breaks), guide="colorbar", na.value="white",name="absorptivity", limits=range(a.breaks))
  
  #----------
  if(scen.k==1) {lper1.1=lper1.map; aper1.1=aper1.map;lper2.1=lper2.map; aper2.1=aper2.map; lper3.1=lper3.map; aper3.1=aper3.map; lper4.1=lper4.map; aper4.1=aper4.map}
  if(scen.k==2) {lper1.2=lper1.map; aper1.2=aper1.map;lper2.2=lper2.map; aper2.2=aper2.map; lper3.2=lper3.map; aper3.2=aper3.map; lper4.2=lper4.map; aper4.2=aper4.map}
  if(scen.k==3) {lper1.3=lper1.map; aper1.3=aper1.map;lper2.3=lper2.map; aper2.3=aper2.map; lper3.3=lper3.map; aper3.3=aper3.map; lper4.3=lper4.map; aper4.3=aper4.map}
  if(scen.k==4) {lper1.4=lper1.map; aper1.4=aper1.map;lper2.4=lper2.map; aper2.4=aper2.map; lper3.4=lper3.map; aper3.4=aper3.map; lper4.4=lper4.map; aper4.4=aper4.map}
  if(scen.k==5) {lper1.5=lper1.map; aper1.5=aper1.map;lper2.5=lper2.map; aper2.5=aper2.map; lper3.5=lper3.map; aper3.5=aper3.map; lper4.5=lper4.map; aper4.5=aper4.map}
  
} #end scen loop

#----------------------------

#FIG. 4 MAPS

#ABS
f1 = plot_grid(aper3.3,aper3.4,aper3.5, ncol = 3, align="v")

abs.plot= grid.arrange(
  arrangeGrob(f1,  
              left=grid::textGrob(label= "latitude (°)", rot=90, gp= gpar(fontsize=18, col="black"))),
  a.leg, 
  widths=c(9,1))

#LAMBDA
f1 = plot_grid(lper3.3,lper3.4,lper3.5, ncol = 3, align="v")

lambda.plot= grid.arrange(
  arrangeGrob(f1,  
              bottom=grid::textGrob(label= "longitude (°)", gp= gpar(fontsize=18, col="black")),
              left=grid::textGrob(label= "latitude (°)", rot=90, gp= gpar(fontsize=18, col="black"))),
  l.leg, 
  widths=c(9,1))
#-------------------

setwd(paste(fdir,"figures/",sep="") )
pdf("Fig4.pdf", height = 6, width = 12)
grid.draw(grobTree(rectGrob(gp=gpar(fill="white", lwd=0)), 
                   grid.arrange(abs.plot, lambda.plot, ncol=1, heights=c(1,1.1) )))
dev.off()  

#********************************************  

#FITNESS OVER TIME WITHOUT EVOLUTION

#check data
sums=colSums(!is.na(lambda.diff[,,1]))
inds= which(colSums(!is.na(lambda.diff[,,1]))>100)
drop.inds= which(colSums(!is.na(lambda.diff[,,1]))<100)

#*******************************************

#calculate fitness differences by scenario
#no plast or evol
lambda.diff= lambda.mean[,,,1 ]

#separate generations
lgen1= lambda.diff[,,1]
lgen2= lambda.diff[,,2]
lgen3= lambda.diff[,,3]

#lambda breaks
lambda.breaks=quantile(lambda.diff, probs = seq(0, 1, length.out=7), na.rm=TRUE )

#scen.labs=  c("plast only", "evol only", "plast + evol", "+ evol of plast")
scen.labs=  c("plast only")

#--------------
#flatten array
l1= array2df(lgen1)
colnames(l1)=c("lambda","year","site")
l2= array2df(lgen2)
colnames(l2)=c("lambda","year","site")
l3= array2df(lgen3)
colnames(l3)=c("lambda","year","site")

#add elevation
l1$elev= pts.sel$elev[l1$site]
l2$elev= pts.sel$elev[l2$site]
l3$elev= pts.sel$elev[l3$site]

#add year
l1$year= years[l1$year]
l2$year= years[l2$year]
l3$year= years[l3$year]

#-------------------
#Figure 5
#surface version of figure

#LAMBDA
for(gen.k in 1:2){
scen.k=1
  
  if(gen.k==1) fcl=l1
  if(gen.k==2) fcl=l2
  
   #drop outliers
    fc= fcl[which(fcl$lambda> (-1) & fcl$lambda< 10),]
    
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
      geom_tile() + theme(legend.position="right")+ 
     scale_fill_distiller(palette="Spectral", na.value="white", name="lambda")+
      theme_classic(base_size=16)+xlab("year")+ylab("Elevation (m)")+
  #  scale_fill_gradientn(colours=brewer.pal(7, "YlOrRd"), values=rescale(lambda.breaks), guide="colorbar", na.value="white",name="lambda")
    
    if(gen.k==1) plot.lambda1= plot.lambda
    if(gen.k==2) plot.lambda2= plot.lambda
    
} #end gen loop
    
    #=============================
 
    #Elevation scatter plots
    
    for(gen.k in 1:2){
    
    #LAMBDA
    #average lambdas across time period
    #first generation and first scenario
    lambda.diff= lambda.mean[,,gen.k,1]
    
    #change erroneous values
    lambda.diff[which(lambda.diff< (-10))]=NA
    lambda.diff[which(lambda.diff> 10)]=NA
    
    per1= colMeans(lambda.diff[which(years %in% 1951:1980),], na.rm = TRUE, dims = 1)
    per2= colMeans(lambda.diff[which(years %in% 1981:2010),], na.rm = TRUE, dims = 1)
    per3= colMeans(lambda.diff[which(years %in% 2011:2040),], na.rm = TRUE, dims = 1)
    per4= colMeans(lambda.diff[which(years %in% 2041:2070),], na.rm = TRUE, dims = 1)
    
    #Set up data
    lper1= cbind(pts.sel, per1)
    lper2= cbind(pts.sel, per2)
    lper3= cbind(pts.sel, per3)
    lper4= cbind(pts.sel, per4)
    
    names(lper1)[9]="lambda"
    names(lper2)[9]="lambda"
    names(lper3)[9]="lambda"
    names(lper4)[9]="lambda"
    
    #elevation plots
    lper1$time.per= "1950-1980"
    lper2$time.per= "1981-2010"
    lper3$time.per= "2011-2040"
    lper4$time.per= "2041-2070"
    
    #plot
    lelev= rbind(lper1,lper2, lper3, lper4)
    
    l.elev= ggplot(lelev) + aes(x = elev, y = lambda)+geom_point()+facet_wrap(~time.per,nrow=2)+geom_smooth(color="gray")+ geom_vline(xintercept=2400, color="gray")+xlab("Elevation (m)")+ylab("Lambda")

    if(gen.k==1)l.elev1=l.elev
    if(gen.k==2)l.elev2=l.elev
    
    } #end gen.k loop
    #-------------------
    
    setwd(paste(fdir,"figures/",sep="") )
    pdf("FigLambda_Scen1.pdf", height = 10, width = 12)
    
    #combine
    plot_grid(l.elev1, plot.lambda1, l.elev2, plot.lambda2, labels = c("A generation 1","B generation 1","C generation 2","D generation 2"), nrow=2, label_size=12, rel_widths = c(1,0.8))
    dev.off()  
    
    #----------------------
    
    