# PLOTS FOR FITNESS SURFACES
#ATTEMPTS NOT USED

# test.spline <- Tps(data.frame(phen$year[site.ind],phen$elev[site.ind]), phen$Jadult[site.ind] )
#  new.grid <- predictSurface(test.spline, nx = 200, ny = 200)
#  image(new.grid)

#------------------
z1= pup.temps["Jadult",inds,,1]
y01= seq(min(pts.sel$elev), max(pts.sel$elev), length.out=length(years) )
fld1<- bicubic(x = years, y = pts.sel$elev[site.ind], z=z1[,site.ind], x0=years, y0=y01 )
gdat <- interp2xyz(fld1, data.frame=TRUE)

plot.Jad= ggplot(gdat) + 
  aes(x = x, y = y, z = z, fill = z) + 
  geom_tile() + 
  geom_contour(color = "white", alpha = 0.5) + 
  scale_fill_distiller(palette="Spectral", na.value="white", name="Jadult") +
  theme_bw(base_size=18)+xlab("year")+ylab("elevation (m)")+theme(legend.position="bottom")

#------------------
#level plots
library(lattice)
levelplot(Jadult ~ elev*year, data = phen, col.regions=terrain.colors(100))

contourplot(Jadult ~ elev*year, data = phen,
            region=TRUE,
            aspect = "fill",
            col.regions = terrain.colors(100))

plot(phen$elev, phen$year)

#subsample
