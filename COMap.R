

##Installs required R packages.
#install.packages("colorRamps")                     

##Loads required R packages

library(raster)
library(maptools)
#library(spgrass6)
library(spatstat)
library(rgdal)
library(RgoogleMaps)
library(colorRamps)
library(dismo) #for extracting google map

library(ggmap)
#_________________________________________________________
#Read points
setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")
pts.sel= read.csv("COpoints.csv")

#_________________________________________________________
#Read lambdas
years= 1950:2099

setwd("C:\\Users\\Buckley\\Google Drive\\Buckley\\Work\\ColiasBiogeog\\OUT\\")
Lambda <- readRDS("lambda.rds")

#Assemble lambdas 
yr.k=1
abs.k=3

#Plot times series
plot(1:50, Lambda[1:50, 2, 3, 1, 1])
#---------------

pts.sel$lambda.gen1= Lambda[yr.k, , abs.k, 1, 1]
pts.sel$lambda.gen2= Lambda[yr.k, , abs.k, 2, 1]
pts.sel$lambda.gen3= Lambda[yr.k, , abs.k, 3, 1]

pts.sel$fat.gen1= Lambda[yr.k, , abs.k, 1, 2]
pts.sel$fat.gen2= Lambda[yr.k, , abs.k, 2, 2]
pts.sel$fat.gen3= Lambda[yr.k, , abs.k, 3, 2]

pts= pts.sel

#-----------------------------------
#set bounding box
dat=pts
dat$lambda= dat$lambda.gen1

#subset to lambda>1
dat2= subset(dat, dat$lambda>1)

#ggmap
bbox <- ggmap::make_bbox(lon, lat, dat2, f = 0.1)
map_loc <- get_map(location = bbox, source = 'google', maptype = 'terrain')

map1 <- ggmap(map_loc, extent='device', base_layer=ggplot(dat2, aes(x=lon, y=lat)))
print(map1)

co.map<- map1 + geom_raster(aes(fill = lambda), alpha=0.5)+ coord_cartesian()
print(co.map)
  #get lat lon back?
