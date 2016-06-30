#DEGREE DAYS CALCULATION

#Single sine wave approximation from Baskerville & Emin 1969
#Double Sine wave approximation of degree days from Allen 1976 
#(see http://www.ipm.ucdavis.edu/WEATHER/ddss_tbl.html)

#LDT is lower developmental threshold

#Calcualted degree days using single sine wave approximation

degree.days=function(Tmin,Tmax,LDT){

dd=NA
if(!is.na(Tmin) & !is.na(Tmax)){

# entirely above LDT
if(Tmin>=LDT) {dd=(Tmax+Tmin)/2-LDT}

# intercepted by LDT
## for single sine wave approximation
if(Tmin<LDT && Tmax>LDT){
alpha=(Tmax-Tmin)/2
#theta1=asin(((LDT-(Tmax+Tmin)/2)/alpha)*pi/180)
#dd=1/pi*(((Tmax+Tmin)/2-LDT)*(pi/2-theta1)+alpha*cos(theta1*pi/180))
theta1=asin(((LDT-(Tmax+Tmin)/2)/alpha))
dd=1/pi*(((Tmax+Tmin)/2-LDT)*(pi/2-theta1)+alpha*cos(theta1))
if(!is.na(dd))if(dd<0){dd=0}
} #matches online calculation

# entirely below LDT
if(Tmax<=LDT){ dd=0}
} #end check NA
return(dd)
}

#--------------------------------------
#RUNS FUNCTION ACROSS MATRIX OF Tmin and Tmax

degree.days.mat=function(Tmat,LDT){

Tmin=Tmat[1]
Tmax=Tmat[2]

dd=NA
if(!is.na(Tmin) & !is.na(Tmax)){

# entirely above LDT
if(Tmin>=LDT) {dd=(Tmax+Tmin)/2-LDT}

# intercepted by LDT
## for single sine wave approximation
if(Tmin<LDT && Tmax>LDT){
alpha=(Tmax-Tmin)/2
#theta1=asin(((LDT-(Tmax+Tmin)/2)/alpha)*pi/180)
#dd=1/pi*(((Tmax+Tmin)/2-LDT)*(pi/2-theta1)+alpha*cos(theta1*pi/180))
theta1=asin(((LDT-(Tmax+Tmin)/2)/alpha))
dd=1/pi*(((Tmax+Tmin)/2-LDT)*(pi/2-theta1)+alpha*cos(theta1))
if(!is.na(dd))if(dd<0){dd=0}
} #matches online calculation

# entirely below LDT
if(Tmax<=LDT){ dd=0}
} #end check NA
return(dd)
}