#Temat is Ta, Rad direct, Rad diffuse, psi

biophys.var_sh=function(Ta, Tg, Tg_sh, u, H_sdir, H_sdif, psi, D, delta, alpha){
  # Ta is ambient temperature, C
  # H_sdir is direct solar radiation flux, W/m^2
  # H_sdif is diffuse solar radiation flux, W/m^2
  # Kt is the clearness index (Wong and chow 2001, Applied energy 69:191-224
  #spec is the row of species data to use
  #SpecDat is the matrix of species data
  #psi is zenith angle, degrees
  
  TaK= Ta+273 #ambient temperature in K
  TaK_sh=TaK
  Tg= Tg+273 #C, T_g- ground surface temperature
  Tg_sh= Tg_sh+273 #C, T_g- ground surface temperature
  u= u *100;  #u- wind speed, convert m/s to cm/s
  H_sdir=H_sdir/10 #divide by ten to convert W/m2 to W/cm2
  H_sdif=H_sdif/10 #divide by ten to convert W/m2 to W/cm2
  
  H_sttl= H_sdir + H_sdif
  
  #Butterfly Parameters
  delta<- delta/10     #delta- thoracic fur thickness, cm
  
  # Biophysical parameters
  r_g=0.30; #substrate solar reflectivity, Kingsolver 1983
  
  #Calculate total surface area as area of cylinder without ends
  A_sttl= pi*D*2 #2 in length  #cm^2
  
  #---------------------------------------------
  
  #Areas, cm^2
  #For basking
  ##A_s,dir, A_s,ref, A_s,ttl- direct, reflected, and total solar radiative heat transfer surface areas 
  A_sdir= A_sttl/2
  A_sref=A_sdir;
  
  #RADIATIVE HEAT FLUx, mW
  Q_s= alpha*A_sdir*H_sdir/cos(psi*pi/180)+alpha*A_sref*H_sdif+alpha*r_g*A_sref*H_sttl;   
  
  #---------------------------------------------		 
  #THERMAL RADIATIVE FLUX
  epsilon_s=0.97; #surface emisivity, ranges from 0.95-1
  sigma= 5.67*10^-9; #Stefan-Boltzman constant, mW cm^-2 K^04 or 5.67*10^-8 Watts m-2 K-4
  
  #Tsky=0.0552*(TaK)^1.5; #Kelvin, black body sky temperature from Swinbank (1963), 
  Tsky= (1.22*Ta -20.4)+273 #K, Gates 1980 Biophysical ecology based on Swnback 1960, Kingsolver 1983 estimates using Brunt equation
  
  Ep=1; #Ep- butterfly thermal emissivity
  #Q_t= 0.5* A_sttl * Ep * sigma * (Tb^4 - Tsky^4) +0.5* A_sttl * Ep * sigma * (Tb^4 - Tg^4)
  
  #---------------------------------------------   	               
  # CONVECTIVE HEAT FLUX
  k_e= 1.3; #k_e- thermal conductivity of the fur, 1.3mWcm^-1*K^-1
  r_i=0.15; #r_i- body radius #Kingsolver 1983
  k_a=0.25; #approximate thermal conductivity of air, mWcm^-1*K^-1, unit conversion checked
  
  v=15.68*10^-2  #cm^2/s, kinematic viscocity of air,  at 300K http://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
  R_e=u*D/v
  N_u=0.6*R_e^0.5
  #N_u=2.3; #Kingsolver 1983;
  
  h_c=N_u*k_a/D;
  h_T=(1/h_c+(r_i+delta)*log((r_i+delta)/r_i)/k_e)^-1;  # h_T- total convective heat tranfer coefficient
  #A_c=A_sttl; #A_c- convective heat transfer surface area
  #Q_c= h_T* A_c* (Tb-Ta);     
  #---------------------------------------------   	 
  #HEAT BUDGET              
  
  # Kingsolver 1983
  #SUN
  #Q_s- total radiative heat flux; Q_t- thermal radiative heat flux; Q_c- convective heat flux
  #Q_s=Q_t + Q_c;
  
  #t solved in wolfram alpha #Solve[a t^4 +b t -d, t]
  a<- A_sttl * Ep *sigma
  b<-h_T * A_sttl
  d<- h_T*A_sttl*TaK +0.5*A_sttl * Ep *sigma*Tsky^4 +0.5*A_sttl * Ep *sigma*(Tg)^4 +Q_s
  
  {Te=1/2*sqrt((2*b)/(a*sqrt((sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)-(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)))-(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)+(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3))-1/2*sqrt((sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)-(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)) }
  
  #SHADE
  #Caclulate without basking by dividing areas by two
  A_sttl=A_sttl/2
  #RADIATIVE HEAT FLUX IN SHADE, mW
  A_sdir= A_sttl/2
  A_sref=A_sdir; 
  H_sdir_sh= 0; #No direct radiation
  H_sdif_sh= H_sdif
  H_sttl= H_sdif + H_sdif_sh #only diffuse and reflected
  Q_s= alpha*A_sdir*H_sdir_sh/cos(psi*pi/180)+alpha*A_sref*H_sdif_sh+alpha*r_g*A_sref*H_sttl; 
  
  #t solved in wolfram alpha #Solve[a t^4 +b t -d, t]
  a<- A_sttl * Ep *sigma
  b<-h_T * A_sttl
  d<- h_T*A_sttl*TaK_sh +0.5*A_sttl * Ep *sigma*Tsky^4 +0.5*A_sttl * Ep *sigma*(Tg_sh)^4 +Q_s
  
  {Te_sh=1/2*sqrt((2*b)/(a*sqrt((sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)-(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)))-(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)+(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3))-1/2*sqrt((sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)-(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)) }
  
  Te=Te-273
  Te_sh= Te_sh-273
  
  #CHECK DATA
  if(!is.na(Te) & !is.na(Te_sh)){
    #select Tb closest to preferred of 35C temperature
    bounds=sort(c(Te, Te_sh))
    all.temps= seq(bounds[1], bounds[2], 0.1) 
    Te= all.temps[which.min(abs(all.temps-(35+273) ))]
  } #end check data
  
  return(Te)
} 

#---------------------------------------------
biophys.var_sh.mat=function(Temat, D, delta, alpha){
  # Ta is ambient temperature, C
  # H_sdir is direct solar radiation flux, W/m^2
  # H_sdif is diffuse solar radiation flux, W/m^2
  # Kt is the clearness index (Wong and chow 2001, Applied energy 69:191-224
  #spec is the row of species data to use
  #SpecDat is the matrix of species data
  #psi is zenith angle, degrees
  Temat= as.numeric(Temat)
  Ta= Temat[1]
  TaK= Ta+273 #ambient temperature in K
  TaK_sh= TaK
  Tg= Temat[2]+273 #C, T_g- ground surface temperature
  Tg_sh= Temat[3]+273 #C, T_g- ground surface temperature
  u= Temat[4] *100;  #u- wind speed, convert m/s to cm/s
  H_sdir=Temat[5]/10 #divide by ten to convert W/m2 to W/cm2
  H_sdif=Temat[6]/10 #divide by ten to convert W/m2 to W/cm2
  psi=Temat[7]
  
  H_sttl= H_sdir + H_sdif
  
  #Butterfly Parameters
  delta<- delta/10     #delta- thoracic fur thickness, cm
  
  # Biophysical parameters
  r_g=0.30; #substrate solar reflectivity, Kingsolver 1983
  
  #Calculate total surface area as area of cylinder without ends
  A_sttl= pi*D*2 #2 in length  #cm^2
  
  #---------------------------------------------
  
  #Areas, cm^2
  #For basking
  ##A_s,dir, A_s,ref, A_s,ttl- direct, reflected, and total solar radiative heat transfer surface areas 
  A_sdir= A_sttl/2
  A_sref=A_sdir;
  
  #RADIATIVE HEAT FLUx, mW
  Q_s= alpha*A_sdir*H_sdir/cos(psi*pi/180)+alpha*A_sref*H_sdif+alpha*r_g*A_sref*H_sttl;   
  
  #---------------------------------------------		 
  #THERMAL RADIATIVE FLUX
  epsilon_s=0.97; #surface emisivity, ranges from 0.95-1
  sigma= 5.67*10^-9; #Stefan-Boltzman constant, mW cm^-2 K^04 or 5.67*10^-8 Watts m-2 K-4
  
  #Tsky=0.0552*(TaK)^1.5; #Kelvin, black body sky temperature from Swinbank (1963), 
  Tsky= (1.22*Ta -20.4)+273 #K, Gates 1980 Biophysical ecology based on Swnback 1960, Kingsolver 1983 estimates using Brunt equation
  
  Ep=1; #Ep- butterfly thermal emissivity
  #Q_t= 0.5* A_sttl * Ep * sigma * (Tb^4 - Tsky^4) +0.5* A_sttl * Ep * sigma * (Tb^4 - Tg^4)
  
  #---------------------------------------------   	               
  # CONVECTIVE HEAT FLUX
  k_e= 1.3; #k_e- thermal conductivity of the fur, 1.3mWcm^-1*K^-1
  r_i=0.15; #r_i- body radius #Kingsolver 1983
  k_a=0.25; #approximate thermal conductivity of air, mWcm^-1*K^-1, unit conversion checked
  
  v=15.68*10^-2  #cm^2/s, kinematic viscocity of air,  at 300K http://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
  R_e=u*D/v
  N_u=0.6*R_e^0.5
  #N_u=2.3; #Kingsolver 1983;
  
  h_c=N_u*k_a/D;
  h_T=(1/h_c+(r_i+delta)*log((r_i+delta)/r_i)/k_e)^-1;  # h_T- total convective heat tranfer coefficient
  #A_c=A_sttl; #A_c- convective heat transfer surface area
  #Q_c= h_T* A_c* (Tb-Ta);     
  #---------------------------------------------   	 
  #HEAT BUDGET              
  
  # Kingsolver 1983
  #SUN
  #Q_s- total radiative heat flux; Q_t- thermal radiative heat flux; Q_c- convective heat flux
  #Q_s=Q_t + Q_c;
  
  #t solved in wolfram alpha #Solve[a t^4 +b t -d, t]
  a<- A_sttl * Ep *sigma
  b<-h_T * A_sttl
  d<- h_T*A_sttl*TaK +0.5*A_sttl * Ep *sigma*Tsky^4 +0.5*A_sttl * Ep *sigma*(Tg)^4 +Q_s
  
  {Te=1/2*sqrt((2*b)/(a*sqrt((sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)-(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)))-(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)+(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3))-1/2*sqrt((sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)-(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)) }
  
  #SHADE
  #Caclulate without basking by dividing areas by two
  A_sttl=A_sttl/2
  #RADIATIVE HEAT FLUX IN SHADE, mW
  A_sdir= A_sttl/2
  A_sref=A_sdir; 
  H_sdir_sh= 0; #No direct radiation
  H_sdif_sh= H_sdif
  H_sttl= H_sdif + H_sdif_sh #only diffuse and reflected
  Q_s= alpha*A_sdir*H_sdir_sh/cos(psi*pi/180)+alpha*A_sref*H_sdif_sh+alpha*r_g*A_sref*H_sttl; 
  
  #t solved in wolfram alpha #Solve[a t^4 +b t -d, t]
  a<- A_sttl * Ep *sigma
  b<-h_T * A_sttl
  d<- h_T*A_sttl*TaK_sh +0.5*A_sttl * Ep *sigma*Tsky^4 +0.5*A_sttl * Ep *sigma*(Tg_sh)^4 +Q_s
  
  {Te_sh=1/2*sqrt((2*b)/(a*sqrt((sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)-(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)))-(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)+(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3))-1/2*sqrt((sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)/(2^(1/3)*3^(2/3)*a)-(4*(2/3)^(1/3)*d)/(sqrt(3)*sqrt(256*a^3*d^3+27*a^2*b^4)+9*a*b^2)^(1/3)) }
  
  Te=Te-273
  Te_sh= Te_sh-273
  
  #CHECK DATA
  if(!is.na(Te) & !is.na(Te_sh)){
    #select Tb closest to preferred of 35C temperature
    bounds=sort(c(Te, Te_sh))
    all.temps= seq(bounds[1], bounds[2], 0.1) 
    Te= all.temps[which.min(abs(all.temps-(35+273) ))]
  } #end check data
  
  return(Te)
} 

#============================

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

#----------------------------
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

#---------------------------

#USE CAMPBELL AND NORMAL MODEL TO ESTIMATE RADIATION
# Biophysical models from Campbell & Norman 1998
# constants

#psi zenith angle radians
#Elevation (m)
#J is Julian Day

calc.rad=function(J.mat, lat=39.74, lon=-105.18, elev=3000){
  
  rho_S=0.7 #rho_S: albedo percent
  
  J= J.mat[1]
  psi= J.mat[2] #zenith angle in radians
  tau= J.mat[3] #transmissivity
  
  sigma=5.67*10^-8 # stefan-boltzman constant, W m^-2 K^-4
  c_p=29.3 # specific heat of air, J/mol degrees K or C
  S_p0=1360 # extraterrestrial flux density, W/m^2 (p159)
  
  # Radiation
  # adjust for elevation
  p_a=101.3* exp (-elev/8200)  # atmospheric pressure
  m_a=p_a/(101.3*cos(psi))  # (11.12) optical air mass
  m_a[which(psi>(80*pi/180))]=5.66
  
  # Flux densities
  #dd2= 1+2*0.1675*cos(2*pi*J/365) #Sears and Angilletta 2012 #dd is correction factor accounting for orbit
  
  #S_p is direct radiation reaching earth's surface
  #S_p=S_p0*tau^m_a*dd2 *cos(psi)  #Sears and Angilletta 2012 #dd is correction factor accounting for orbit
  S_p=S_p0*tau^m_a *cos(psi) #Use initial Cambell and Norman
  
  #S_d=0.3*(1-tau^m_a)* S_p0*cos(psi)*dd2 #with correction factor
  S_d=0.3*(1-tau^m_a)* S_p0*cos(psi)
  #S_t=S_p*cos (psi)+S_d # solar irradience 
  S_r= rho_S*(S_p+S_d) # (11.10) reflected radiation
  
  #return direct, diffuse, reflected
  return( c(S_p, S_d, S_r))
}

#-------------------------

#https://www.r-bloggers.com/approximate-sunrise-and-sunset-times/
# OR suncalc getSunlightTimes

suncalc<-function(d,Lat=48.1442,Long=-122.7551){
  ## d is the day of year
  ## Lat is latitude in decimal degrees
  ## Long is longitude in decimal degrees (negative == West)
  
  ##This method is copied from:
  ##Teets, D.A. 2003. Predicting sunrise and sunset times.
  ##  The College Mathematics Journal 34(4):317-321.
  
  ## At the default location the estimates of sunrise and sunset are within
  ## seven minutes of the correct times (http://aa.usno.navy.mil/data/docs/RS_OneYear.php)
  ## with a mean of 2.4 minutes error.
  
  ## Function to convert degrees to radians
  rad<-function(x)pi*x/180
  
  ##Radius of the earth (km)
  R=6378
  
  ##Radians between the xy-plane and the ecliptic plane
  epsilon=rad(23.45)
  
  
  ##Convert observer's latitude to radians
  L=rad(Lat)
  
  ## Calculate offset of sunrise based on longitude (min)
  ## If Long is negative, then the mod represents degrees West of
  ## a standard time meridian, so timing of sunrise and sunset should
  ## be made later.
  timezone = -4*(abs(Long)%%15)*sign(Long)
  
  ## The earth's mean distance from the sun (km)
  r = 149598000
  
  theta = 2*pi/365.25*(d-80)
  
  z.s = r*sin(theta)*sin(epsilon)
  r.p = sqrt(r^2-z.s^2)
  
  t0 = 1440/(2*pi)*acos((R-z.s*sin(L))/(r.p*cos(L)))
  
  ##a kludge adjustment for the radius of the sun
  that = t0+5 
  
  ## Adjust "noon" for the fact that the earth's orbit is not circular:
  n = 720-10*sin(4*pi*(d-80)/365.25)+8*sin(2*pi*d/365.25)
  
  ## now sunrise and sunset are:
  sunrise = (n-that+timezone)/60
  sunset = (n+that+timezone)/60
  
  return(list("sunrise" = sunrise,"sunset" = sunset))
}
