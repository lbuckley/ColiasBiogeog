#Estimate air pressure in kPa #http://www.engineeringtoolbox.com/air-altitude-pressure-d_462.html
airpressure_elev<- function(h){  #H is height in meters 
  p= 101325* (1 - 2.25577*10^(-5)*h)^5.25588       
  p= p/1000 #convert to kPa
  return(p)
}

# V_r is wind velocity at reference height
V_z<- function(V_r, z_0=0.02, z_r=1.54, z=0.2){ V_r*log((z+z_0)/z_0+1)/log((z_r+z_0)/z_0+1)}
#from Porter et al 1973

#z_0 is surface roughness
#z_r is initial reference height
#z is height to scale to

#SEE SurfaceRoughness_12May2014.R for C1 surface roughness estimate
#z_0=0.02 #m

#calculate air temp at some height z
air_temp_at_height_z<-function(z_0, z_r, z, T_r, T_s){
  T_r= as.numeric(T_r)
  T_s= as.numeric(T_s)
  T_z<-(T_r-T_s)*log((z+z_0)/z_0+1)/log((z_r+z_0)/z_0+1)+T_s ##this is exactly eqn (19) of the notes
  return(T_z)
}

air_temp_at_height_z_mat<-function(Tdat, z_0, z_r, z){
  T_r=Tdat[1]; T_s=Tdat[2]
  T_z<-(T_r-T_s)*log((z+z_0)/z_0+1)/log((z_r+z_0)/z_0+1)+T_s ##this is exactly eqn (19) of the notes
  return(T_z)
}

#---------------------------------------

### TEMP ACROSS HOURS
Thours.mat=function(Tmat, alpha=1.86, beta= -0.17, gamma=2.20){
  #Tmx= max temperature
  #Tmn= min temperature
  
  T=NA
  if( sum(is.na(Tmat))==0){
    
    Tmx= as.numeric(Tmat[1])
    Tmn= as.numeric(Tmat[2])
    tr= as.numeric(Tmat[3])
    ts= as.numeric(Tmat[4])
    
    l= ts-tr #daylength
    
    tx= 0.5*(tr+ts)+alpha #time of maximum temperature
    tn= tr+ beta #time of minimum temperature
    
    hrs= 1:24
    
    Tsn= Tmn+(Tmx-Tmn)*sin((pi*(ts-tr-beta))/(l+2*(alpha-beta)))
    
    #daytime
    T= Tmn+(Tmx-Tmn)*sin((pi*(hrs-tr-beta))/(l+2*(alpha-beta)))
    #am
    inds= which(hrs<= (tr+beta))
    T[inds]= Tmn+(Tsn-Tmn)*exp(-(gamma*(hrs[inds]+24-ts))/(24-l+beta))
    #pm
    inds= which(hrs>ts)
    T[inds]= Tmn+(Tsn-Tmn)*exp(-(gamma*(hrs[inds]-ts))/(24-l+beta))
    
  } #end check NA
  return(T)
}

#-----------------------
#RADIATION WITH CLOUDINESS, returns direct and diffuse radiation
library(ks)

Rad.mat=function(rad){
  #PICK TAU, atmospheric transmisivity, from distribution for k_t
  #USES KERNAL ESTIMATE FIT TO HOURLY NREL SRRL DATA (loaded when sourcing solar radiation functions)

  pick.tau= function(hr) rkde(fhat= kdes[[ round(hr)-5]], n=1 )
  
  hrs=1:24
  taus= rep(0, 24)
  taus[6:20]= apply(matrix(6:20, ncol=1),FUN=pick.tau, MARGIN=1)
  
  #constraint
  taus[taus<0]=0
  #CONSTRAIN TAUs BETWEEN 7AM AND 3PM,
  seq1= 7:15
  taus[seq1[which(taus[7:15]<0.2)]]=0.2
  
  #USE ERBS TO REPARTITION RADIATION (Olyphant 1984)
  #Separate Total radiation into components
  #kt is clearness index
  # Models presented in Wong and Chow. 2001. Applied Energy 69(2001):1991-224
  #Use Erbs et al model
  
  #kd- diffuse fraction
  kd= rep(NA, length(taus))
  
  inds= which(taus<0.22) 
  kd[inds]= 1-0.09*taus[inds]
  inds= which(taus>0.22 & taus<0.8) 
  kd[inds]= 0.9511 -0.1604*taus[inds] +4.388*taus[inds]^2 -16.638*taus[inds]^3 +12.336*taus[inds]^4
  inds= which(taus>=0.8)
  kd[inds]= 0.125 #Correction from 16.5 for Niwot from Olyphant 1984
  
  #return direct and diffuse
  rad= as.numeric(rad)
 return (c(rad*(1-kd),rad*(kd)) )
}  

#----------------------------
LambdaCalc= function(Eggs, SurvMat, SurvDaily, MaxEggs)

if(!is.nan(Eggs)){
  MaxDay=5
  Lambda1=0
  for(day in 1:MaxDay){
    Eggs1= min(Eggs, MaxEggs-Eggs*(day-1))  ###LIMIT MAX NUMBER EGGS
    if(Eggs1<0) Eggs1=0
    Lambda1= Lambda1+ SurvMat * SurvDaily^day *Eggs1;                        
  }#end loop days
 return(Lambda1) 
}