#DDD model Module based

#----------------------------------------------------------------------------------
#     Script name:  DDDModulFunc
#     Author:     Thomas Skaugen & Zelalem Mengistu
#     Revised: 17.12.2019
#              
#Purpose   : DDD (Distance Distribution Dynamics) hydrological model simulates
#               1)saturated soil water flow, states of saturated (2-D)) and unsaturates zones
#               2) snow accumulation, distribution and melt
#               3) river run off for a given catchment.
#               4)saves the states for each day where precip > 0, i. e where SCHADEX might insert simulated extreme precip
#
# Model input: 1) catchment physical/geographical description
#              2) weather data, temperature,precipitation and river discharge
#              3) model states variables (optional)
#
# Model output: simulated 1)saturated and unsaturated soil water
#                         2)snow water equivalent and snow distribution
#                         3)river discharge
#                         4)state files for each rainy day. Remember to include temperature for 5 days t-1, t=0, t+1, t+2, t+3 
#
# Running DDD:    It is possible to save model state variables and run the model
#                 starting from saved state variables
#----------------------------------------------------------------------------------


function DDDModulFunc(Gpar, startsim, tprm, prm, ptqfile, utfile, r2fil, modstate, savestate,catchment, kal)
                                 

#-------------------------------------preprocess start------------------------------------------
#-----------------------------------------------------------------------------------------------

#####################VARIABLES FOR 10 ELEVATION ZONES###########################################
isoil = zeros(10)              # precipitation and snowmelt from the elevation zones
gwgt = zeros(10)               # weights for glaciermelt pr each elevation zone
swgt = zeros(10)               # weights for input to soils pr each elevation zone NB bogs are not part of this yet
gisoil = zeros(10)             # glaciermelt from the elevation zones
spd = zeros(10)
swe_h = zeros(10)        #SWE pr altitude level
wcd = zeros(10)
sca = zeros(10)
nsno = zeros(10)
alfa = zeros(10)
ny = zeros(10)
hprecip = zeros(10)
htemp = zeros(10)
hfelt = zeros(10)
snowfree = zeros(10) #Indicator variable 1 is sca  less than gca 0 if sca > gca
tempvec = zeros(10) # temperature vector  to save in statefile
PR = zeros(10)
PS = zeros(10)
MWEB = zeros(10)
MWCX = zeros(10)
MW = zeros(10)
MWGLAC = zeros(10)
MWGLAC1 = zeros(10)    
scaobx = zeros(10)
SWrad = zeros(10)
LA = zeros(10)
LT = zeros(10)
SH = zeros(10)
LE = zeros(10)
A = zeros(10)
Pa = zeros(10)
Cl = zeros(10)    

#Addition for Enegy balance snowmelt
snro = zeros(10)            #Snow density
melt = zeros(10)            #instead of degreeday melt (MW)
snowdepth = zeros(10)        #hmm..
sndens = zeros(10)           #hm ...
#initializing

#Reading parameters
hfelt[1] = prm.val[3]+(prm.val[4]-prm.val[3])/2.0
hfelt[2] = prm.val[4]+(prm.val[5]-prm.val[4])/2.0
hfelt[3] = prm.val[5]+(prm.val[6]-prm.val[5])/2.0
hfelt[4] = prm.val[6]+(prm.val[7]-prm.val[6])/2.0
hfelt[5] = prm.val[7]+(prm.val[8]-prm.val[7])/2.0
hfelt[6] = prm.val[8]+(prm.val[9]-prm.val[8])/2.0
hfelt[7] = prm.val[9]+(prm.val[10]-prm.val[9])/2.0
hfelt[8] = prm.val[10]+(prm.val[11]-prm.val[10])/2.0
hfelt[9] = prm.val[11]+(prm.val[12]-prm.val[11])/2.0
hfelt[10] = prm.val[12]+(prm.val[13]-prm.val[12])/2.0
phi = prm.val[15]*pi/180  #converts Lat degrees to radians lat
thi = prm.val[16]*pi/180  #converts Lon degrees to radians Lon 
pkorr = tprm[1]       # prm.val[17]   # pskorr,prkorr   :Met corrections
skorr = tprm[2]       # prm.val[18]   # pskorr,prkorr   :Met corrections
u =  tprm[3]          # prm.val[19]   # wind
pro = tprm[4]         # prm.val[20]   # % liquid water content in snow
TX = tprm[5]          # prm.val[21]   # Threshold temp for rain/snow
CGLAC = prm.val[22]       # prm.val[22]    # CGLAC          :Degreeday index for glaciermelt
a0 = tprm[6]# prm.val[23]  #                 :Alfa null in snofall unit distribution
d = prm.val[24]   #                 :Decorrelation length (Measured in n) for precipitation
CX = tprm[7]               # prm.val[25]  # degreeday factor
TS = tprm[11]              #prm.val[26]    # Threshold temperature for snow and glacial melt
Timeresinsec = prm.val[27] # in seconds
MAD1 = prm.val[28]          # mean annual discharge
totarea = prm.val[29]      # total area in m2
maxLbog = prm.val[30]     #Maximum distance in wetlands
midLbog = prm.val[31]      # mean distance in wetlands
bogfrac = prm.val[32]      #areal fraction wetlands
zsoil = prm.val[33]        #frac zero dist soil
zbog = prm.val[34]         #frac zero dist bog
NoL = Int(prm.val[35])          #Number of layers
R = prm.val[36]            #soil moisture content/100 for field capacity of soils 
CritFlux = prm.val[37] # critical flux for initiating Dynamic River Networks
GshInt = tprm[9]# prm.val[38]       # shapeparameter saturation Capital lambda
GscInt = tprm[10]# prm.val[39]       # scaleparameter saturation Capital lambda
Gshape, Gscale = Big2SmallLambda(GshInt, GscInt) # Coverting integrated celerity to layers    
#Gshape = tprm[9] #Gpar[1]
#Gscale = tprm[10]#Gpar[2]
rv = tprm[8]              # prm.val[40]           # celerity of water in rivernetwork
midFl = prm.val[41]        # mean of distance in river network
stdFl = prm.val[42]        # standard deviation of distances in river network
maxFl = prm.val[43]         # maximum distance in river network
maxDl = prm.val[44]        # maximum distance in hillslopes
midDL = prm.val[45]       # mean distance in hillslopes
glacfrac = prm.val[46]    # fraction of glaciers
midGl = prm.val[47]       # mean of distances of glaciers
stdGl = prm.val[48]       # standard deviation of distances of glaciers
maxGl = prm.val[49]       # maximum distance of glaciers
g1    = prm.val[50]       # areal fraction of glaciers in first elevation zone
g2    = prm.val[51]
g3    = prm.val[52]
g4    = prm.val[53]
g5    = prm.val[54]
g6    = prm.val[55]
g7    = prm.val[56]
g8    = prm.val[57]
g9    = prm.val[58]
g10   = prm.val[59]
Lakefrac = prm.val[60]     # effective lake % 
Lv = prm.val[61]           # lake celrity [m/s] Can it be fixed?
midLake = prm.val[62]      # mean distance to outlet in Lake
stdLake = prm.val[63]      # standard deviation of distances in Lake 
maxLake = prm.val[64]      # maximum distance in Lake
snfjell = 100.0*prm.val[65]# in parameterfile as fraction 
meandailyP = prm.val[66]   # daily mean value
meandailyT = prm.val[67]   # daily mean value


#Merging assumed Glacier river network  with observed river network. We assume independent normal distributions
maxFl = maxFl + maxGl
midFl = midFl + midGl  #mean of flowlength distribution
stdFl = stdFl + stdGl

     # reading ptq data, elevation zones
if (Timeresinsec >= 3600)    
  ptqinn= CSV.read(ptqfile,DataFrame,header=["yr","mnt","day","hr","p01","p02","p03","p04","p05","p06","p07","p08","p09","p10",
                  "t01","t02","t03","t04","t05","t06","t07","t08","t09","t10","q"], delim=';')           
else
  ptqinn = CSV.read(ptqfile,DataFrame,header=["yr","mnt","day","hr","min","p01","p02","p03","p04","p05","p06","p07","p08","p09","p10",
                  "t01","t02","t03","t04","t05","t06","t07","t08","t09","t10","q"],delim=';') 
end        
days = Int(length(ptqinn.yr))    

MAD2 = exp(-4.59 + 1.135*log(meandailyP)+0.97*log(totarea/1000000)-0.0603*meandailyT + 0.053*log(snfjell)-0.00014*hfelt[5])  # R2 = 0.99, Area in km2, Snfjell in fraction   
MAD = MAD1
if (glacfrac > 0.1 ) 
     MAD = MAD2    # needs to adjust MAD if glaciers are present since MAD is used to estimate the groundwater reservoir
end 

#Constants
hson = 10                               # Elevation zones
unitsnow = 0.1                          # mean of unit snow [mm]
n0 = unitsnow*a0                        # shape parameter of unit snow
gtcel = 0.99                            # threshold groundwater capacity -> Overland flow
CFR = 2.5*(Timeresinsec/86400)*0.0833   # Fixed as 1/12 of estimate of CX= 2.5 for 24 hours 
len = Int(5*(86400/Timeresinsec))       # number og timestepes to estimate snowpack temperature. first factor is days)
startsim = startsim + len               # taking into account estimating snowpack temperature
STempvec = zeros(len)                   # Temperature vector for estimating snowpack temperature
taux = 0.0

Pa[1:10] = 101.3*((293 .- 0.0065.*hfelt[1:10])/293).^5.26     # Air pressure as a function of height Zhu et al. 2012 and Stoy et al, 2018
MPa = mean(Pa)

#vectors and matrixes
ddist = zeros(NoL)                       # distribution of water in layers
ddistx = zeros(NoL)                      # water present in layers
ddistxtemp = zeros(NoL)                  # temporay of the above
aktMag = zeros(NoL)                      # current water in Layers
qberegn = zeros(days)                    # for calculation of skill scores etc.
nodaysvector = zeros(Int64,NoL)

if(kal == 0) 
    if (Timeresinsec >= 3600)    
      simresult = zeros(days,35)            # matrix which into the results are written
    else
      simresult = zeros(days,36)
    end
end  

#if(ebelements==1) energyel <- matrix(0.0, (days-startsim+1),194) #matrix to which the simulated energybalance elements are written
#grwpointresult = zeros(days,26)    # matrix which into the groundwater point results are written

# Areas and weights
area = (1-(bogfrac+glacfrac))*totarea    #area not covered by wetlands or glaciers. These three landscape types have different distancedistributions
area2 = (1-(bogfrac))*totarea           #area in which we have hillslope process and rivernetwork processes
bogarea = bogfrac*totarea               #area af wetlands
glacarea = glacfrac *totarea             #area of glaciers
elevarea = (totarea/hson)               #area pr elevationzone
gca = [g1,g2,g3,g4,g5,g6,g7,g8,g9,g10] #fraction of glacier pr elevation zone 
soilca = 1 .- gca

# gwgt is the fraction of glaciated area pr elevation zone in relation to total glaciated area 
if glacarea > 0
    gwgt = gca .* elevarea/glacarea
    else
    gwgt .= 0
end
# Finds the fraction of soils (and bogs) in each elevation zone in relation to total soil (and bog) area 
swgt = soilca .* elevarea/(area + bogarea)

#Subsurface celerities 
k = CeleritySubSurface(NoL, Gshape, Gscale, midDL, Timeresinsec)  # Celerity of subsurface (and overland) flow
    
#Surfacecelerities
k[1]  = 0.007                             # Overland flow permeable areas celerity from Litterature 
bogspeed = k[1]*1.0                       # Same celerity as overland flow

#Unit hydrographs for Wetlands
antBogsteps = Int(trunc(maxLbog/(bogspeed*Timeresinsec))+1) # Number of timsteps to drain wetlands
UHbog = zeros(antBogsteps)                                  # Vector for Bog unit hydrograph
if (bogarea > 0) 
    UHbog = SingleUH(bogspeed,Timeresinsec, midLbog, maxLbog, zbog)
end    
if(bogarea <= 0 )
    UHbog = 1.0
end

#Unit hydrographs for hillslopes based on celerities and distance distributions,exponentially distributed
antHorlag = Int(trunc(maxDl/(k[NoL]*Timeresinsec))+1)       # +1 includes the final day.  NB only in soils part
nodaysvector = zeros(Int16, NoL)                                   # Number of days to drain each layer
nodaysvector[1:NoL] = Int.(trunc.(maxDl./k[1:NoL]./Timeresinsec) .+1)   # integer time steps

layerUH = zeros(NoL,antHorlag)

for i in 1: NoL
    layerUH[i,1:Int(nodaysvector[i])] .= SingleUH(k[i], Timeresinsec, midDL, maxDl, zsoil)
end

#UH for RIVER, normally distributed 
#number of time units in river routing, rounded up. If set equal to one day, time is actually less
UHriver, noDT, nodaysRiv = SingleNormalUH(rv,Timeresinsec,midFl,stdFl,maxFl) #This subroutine extracts a UH normally distributed. 

#UH for LAKES, normally distributed 
#nodaysLake = Int(trunc(maxLake/(Lv*Timeresinsec))+1) # Number of timsteps to drain Lakes Needed if exponential UH
#noDTLake = nodaysLake                                # Needed if exponential UH
#number of time units in Lake routing, rounded up. If set equal to one day, time is actually less
UHLake, noDTLake, nodaysLake = SingleNormalUH(Lv,Timeresinsec,midLake,stdLake,maxLake) #This subroutine extracts a UH normally distributed. 
#UHLake = SingleUH(Lv,Timeresinsec,midLake,maxLake,0) #This subroutine extracts an UH exponentially distributed. 
 
vectorlengde = Int(antHorlag+noDT-1)
qsimutx = zeros(noDT)          # vector for runoff i m3/s
qsimut = zeros(noDT)           # vector for runoff in mm
qBogut = zeros(noDT)           # vector for runoff from wetlands bogs
QRivx = zeros(noDT)            # all qsimutvectors are stored antHorlag timesteps back
QLakex = zeros(noDT)           # all qsimutvectors are stored antHorlag timesteps back

#Groundwater layers; 2dim levels, 1 fastest, NoL slowest#

Layers = zeros(NoL, antHorlag)
BogLayers = zeros(antBogsteps)
LakeLayers = zeros(nodaysLake)
    
Magkap, M = LayerEstimation(GshInt,GscInt,Timeresinsec,maxDl,midDL, MAD, totarea, NoL, gtcel) 
    
#Extracting areas and heights in pyramide-plot
# Areas: Important: same dimensions as Layers. The area drained pr time-step box for each layer
Areas, delta_d = PyrAreas(NoL,totarea,maxDl,nodaysvector, layerUH, antHorlag) # A matrix of areas[in m2] for each time-step box, 

#States
modstate = 0
#if(modstate ==1) load(tilst.file)  # Run with state
if (modstate == 0)                   # Run with initial values  
    sm = 0.0                     # Soilmoisture
    smbog = 0.0                   # Saturation of Bogs
    
    # groundwater initial value
    for i in 1:NoL
        Layers[i,1:nodaysvector[i]] .= 0.0 ./ nodaysvector[i]
     end 
    QRivx[1:noDT] .= qsimutx                             #discharge in m3/s
    totdef = M
end


isoilx = 0.0
Snowdepth = 0.0
snittT = 0.0
ea = 0.0
ea_S = 0.0
                                 
############################################################################
#                                                                          #
#simulation starts here.............. days is no. time steps               #
############################################################################
for i in startsim:days

dato = Date(ptqinn.yr[i], ptqinn.mnt[i] ,ptqinn.day[i])                #date
    
hr =ptqinn.hr[i]  
    
#if (Timeresinsec == 86400)
#   ptqinn.hr[i] = 12 
#end
DN = Dates.dayofyear(dato)                                             #daynumber 
#println("Timestep= ", i," max= ",days, " dato= ",dato," DN= ", DN)  
 
#Reads Precipitation and temperature for each elevation zone  
if (Timeresinsec >= 3600)         
  htemp = ptqinn[i, 15:24]  #<-FromReadPrecipTemp$htemp
  hprecip = ptqinn[i,5:14]   #<- FromReadPrecipTemp$hprecip
else
  htemp = ptqinn[i, 16:25]  #<-FromReadPrecipTemp$htemp
  hprecip = ptqinn[i,6:15]   #<- FromReadPrecipTemp$hprecip
end

meanprecip = mean(hprecip)
meantemp = mean(htemp)

##########Snow simulations#####################################################
for idim in 1:hson               #elevation zones

if (Timeresinsec >= 3600)         
  STempvec = TemperatureVector(ptqinn[(i-len+1):i,15:24],idim, len)
else
  STempvec = TemperatureVector(ptqinn[(i-len+1):i,16:25],idim, len)
end
  
  
#Estmates the mean snowpack temperature
     snittT = SnowpackTemp(STempvec) 
        
#Calculates Rain vs Snow and Snow/Glaciermelt
WaterIn = NedbEB(CFR,DN,hprecip[idim],htemp[idim],spd[idim],Timeresinsec,TX,pkorr,skorr,hr,u,Pa[idim],CGLAC,CX,TS, taux,snittT,phi,thi)
               
#From WaterIn per elevation zone: 
        #PS,PR,MW,MWGLAC,SWrad, LA,LT,SH,LE,GH,PH,CC,A,taux,snittT, Tss,Cl, RH
  PS[idim] = WaterIn[1]         # Precipitation as snow
  PR[idim] = WaterIn[2]         # Precipitation as rain
  MWEB[idim] = WaterIn[3]         # Snow melt potential from Energy balance
  MWGLAC[idim] = WaterIn[4]   # Glacier melt potential
  SWrad[idim] = WaterIn[5]    # Short wave radiation (net)
  LA[idim] = WaterIn[6]        # Long wave radiation atomspheric
  LT[idim] = WaterIn[7]          # Long wave radiation terrestrial
  SH[idim] = WaterIn[8]          # Turbulent heat
  LE[idim] = WaterIn[9]          # Latent heat
  GH = WaterIn[10]          # Ground heat
  PH = WaterIn[11]          # Precipitation heat
  CC = WaterIn[12]          # Cold content
  A[idim] = WaterIn[13]            # Albedo
  taux = WaterIn[14]     # Aging factor of snow for albedo calulations
  snittT = WaterIn[15]  # Snowpack temperature
  Tss = WaterIn[16]        # Snow surface temperature
  Cl[idim] = WaterIn[17]          # CloudCover
  RH = WaterIn[18]         # Relative humidity
  MWGLAC1[idim]=WaterIn[19]	    #Potential glacial melt from EB 
  MWCX[idim] = WaterIn[20]   	#Potential snowmelt from degreeday factor
            
  gisoil[idim] = MWGLAC1[idim]
  
#######################
#  here you choose whether snowmelt is to be carried out by CX (degree day factor) or EB (Energy balance)  
   #MW[idim] = MWCX[idim]  # potential melt as from CX (degreeday) estimate
    MW[idim] = MWEB[idim]  # potential melt as from EB estimate
#######################

# calculates Snow accumulation and Snowdistribution            

FromSnow = SnowGamma(PR[idim],PS[idim],MW[idim],sca[idim],spd[idim],wcd[idim],pro,nsno[idim],alfa[idim],ny[idim],a0,n0,a0,d)
    isoil[idim] = FromSnow[1]  #precipitation and snowmelt
      spd[idim] = FromSnow[2] 
      wcd[idim] = FromSnow[3] 
      sca[idim] = FromSnow[4] 
     nsno[idim] = FromSnow[5] 
     alfa[idim] = FromSnow[6]
       ny[idim] = FromSnow[7]
             
    if (spd[idim] > 10000)     
	 if(kal == 1)
       spd[idim]  = 8000.0
	   skorr= 0.9*skorr
       println("Skorr is reduced, you build snowtowers = trend in SWE")
	  end 
    end
    
    if (sca[idim] < gca[idim])
      snowfree[idim] = 1.0
    else
      snowfree[idim] = 0.0                   #test for glaciermelt,  no melt if snowcovered, snowfree[idim] =0.0  
    end
   
    #Snowdepths and snowdenities
    snowdepthx = snowdepth[idim]
    snomag =  (sca[idim]*spd[idim])           #snowreservir this timestep per idim

    # Snowdepth, Snowdensity 
    snowdepth[idim],sndens[idim] = NewSnowSDEB(PR[idim],PS[idim],htemp[idim],TX,MW[idim],snomag,snomag,snowdepthx)
         
  end # for elevation zones

  MSWrad = mean(SWrad)
  MLA = mean(LA)
  MLT = mean(LT)
  MSH = mean(SH)
  MLE = mean(LE)
  MALB =mean(A)
  midGLAC = mean(MWGLAC)  # degree-day
  midGLAC1 = mean(MWGLAC1) # EB 
  midMWEB = mean(MWEB) # mean snowmelt from EB
  midMWCX = mean(MWCX) # mean snowmelt from CX
  mCl = Cl[5] 

  #Snowreservoir
  snomag = mean(sca[1:hson] .* spd[1:hson]) #mean catchment SWE ,must multiply with SCA to have arealvalues

  swe_h[1:hson] = sca[1:hson] .* spd[1:hson] #SWE pr. elevation zone

  middelsca = sum(sca[1:hson])/hson

  snofritt = 1-middelsca
  wcs = mean(wcd[1:hson])
  
  #Estimate isoil from all elevation zones
  misoil = sum(isoil.*swgt)                                      #snowmelt and rain on soils. 
                                                                 # isoil is weighted by the fraction of soils
                                                                 # in relation to soil an bog area pr elevation band 
  m_r_onglac = sum(isoil.*gwgt)                                  # snowmelt and rain from glaciated area. 
                                                                 # isoil is weighted by the fraction of glaciers 
                                                                 # pr elevation band in relation to total glaciated area  
                                                                 # glacier melt (gisoil) in mm this timestep. 
                                                                 # m_r_onglac + outglac is total output from glacier
  outglac = sum(gisoil[1:hson].*gwgt[1:hson].*snowfree[1:hson])  # gwgt because it is going to be scaled by glacfrac later on

  #Estimating current capacity in Layers
  ddistx = LayerCapacityUpdate(Layers, nodaysvector, Magkap, NoL)
  
  # State of unsaturated zone, it can etiher be zero (complete saturation) or positive, the deficit. If negative we have overland flow
  totdef = sum(ddistx[2:NoL])       # This may be negative if input outstrips deficits
  D = max(totdef,0)

  toSoil = misoil*(1-(glacfrac))+glacfrac*(outglac+m_r_onglac) # input from rain, snow and glaciers always > 0    
 
  PotEvap = PotentialEvapPT(meantemp, MSWrad, MLA, MLT, MPa)
  if(PotEvap < 0.0)
    PotEvap = 0.0
  end

  #calculating evapotranspiration from soilmoisture reservoir
  sm, ea_sm, ea_S = UnsaturatedEvapEB(toSoil, meantemp, sm, M, D, PotEvap)
                  # updated sm
                  # evapotranspiration to be drawn from the soilmoisture
                  # evapotranspiration to be drawn from Layers    

  #calculating additional (ea_S) evapotranspiration directly from Layers and updating Layers
  Layers, ea_S = LayerEvap(Layers, nodaysvector, ea_S, layerUH, NoL)

  ea = ea_sm + ea_S # adding evapotranspiration from sm (first) and Layers (second)
  
  #Checking for capacity after evapotranspiration
  ddistx = LayerCapacityUpdate(Layers, nodaysvector, Magkap, NoL)
 
  # State of unsaturated zone, after Layers are updated (after evapotranspiration) 
  totdef = sum(ddistx[2:NoL])       # This may be negative if input outstrips deficits
  D = max(totdef,0)
 
# Call for the soilwater routine  calculates soil water (unsaturated)
  outx, sm = UnsaturatedExEvap(toSoil, sm, R, D)

  #From Wetlands/Bogs
  FromWetlands = WetlandsEB(misoil, meantemp, middelsca, smbog, M, PotEvap)
  outbog = FromWetlands[1]
   smbog = FromWetlands[2]
   eabog = FromWetlands[3]
 
  #Estimating input to Layers according to capacity
  ddist = GrvInputDistribution(outx,NoL,ddistx, ddist)
 
 OFDDD = 1
  #Calculating UH for overland flow Layer #1
  if (OFDDD == 1)
       
     if(ddist[1]*outx > 0.0)# overland flow is confirmed for this event
        #println("before DRN",layerUH[1,1:nodaysvector[1]])
        layerUH[1,1:nodaysvector[1]] = OverlandFlowDynamicDD(k,ddist,outx, layerUH, nodaysvector,NoL, 
                    midDL, CritFlux, Timeresinsec)
        #println("after DRN",layerUH[1,1:nodaysvector[1]])
     end
 end
     
 GDT = sum(Layers[1:NoL,1]) + sum(ddist .* outx .* layerUH[1:NoL,1]) #Groundwater to be discharged into the river network + this timesteps contribution
 GDT_Bog = BogLayers[1] + outbog*UHbog[1]            #Bogwater to be discharged into the river network + this timesteps contribution
    
  #Updating the saturation Layers
  LayerUpdate!(ddist, outx, Layers, layerUH, nodaysvector, NoL) # mm
  BogLayerUpdate!(outbog, BogLayers, UHbog, antBogsteps)     # mm
   
  lyrs = sum(Layers)# gives the sum of layers after todays runoff has taken place
  boglyrs = sum(BogLayers)
  
  #reinstate original overland flow layer        
  layerUH[1,1:nodaysvector[1]] = SingleUH(k[1], Timeresinsec, midDL, maxDl, zsoil)
        
  #GrWPoint <- GrW_Point(nodaysvector, i, NoL, Areas, Layers, totarea, mindist[RN,],delta_d)
  
  #Response from the soils part [m3/s] from todays input. This is a vector giving todays runoff and antHorlag + nodaysRiv ahead
  qsimutx .= (((GDT/1000)*area2)/Timeresinsec).* UHriver    #from mm to m3/s       
    
  #Response from Wetlands in m3/s
  if(bogarea > 0.0 )
     qBogut .= (((GDT_Bog/1000)*bogarea)/Timeresinsec) .* UHriver #from mm to m3/s
  else
     qBogut[1:noDT] .= 0.0  #m3/s
  end 

  #Total hillsope response
  qsimutx[1:noDT] .= qsimutx[1:noDT] .+ qBogut[1:noDT]     #adding contribution from bogs all in m3/s
  
  #Updating the routing in river network
  QRivx, QRD = RiverUpdate(noDT,QRivx,qsimutx)
  qmm_state = (sum(QRivx) - QRD)*(Timeresinsec*1000/totarea)       # This is also a reservoir [mm], water from todays event is stored for future runoff in the RN, relevant for WB.  
  
  GDT_Lake = LakeLayers[1] + QRD*UHLake[1]            #Lakewater to be discharged from outlet + this timesteps contribution, i.e. catchemnt response 
  BogLayerUpdate!(QRD, LakeLayers, UHLake, nodaysLake)#
  Qm3s = GDT_Lake # runoff in m3/s
  Qmm = (GDT_Lake*Timeresinsec*1000/totarea)  #GDT_Lake in mm/Timeresinsec
  
  fac = 1000.0/Timeresinsec
  if (kal==0)
    if (Timeresinsec >=3600)
  #Assigning outdata to vector, one vector for each timestep
   simresult[i, 1:4] .= [ptqinn.yr[i],ptqinn.mnt[i], ptqinn.day[i], ptqinn.hr[i]]     
   simresult[i,5:8] .= [round(meanprecip,digits=4), round(meantemp,digits=3), ptqinn.q[i], round(Qm3s,digits=6)]
   simresult[i,9:12] .= [round(middelsca,digits=3),round(snomag,digits=4),round((M-totdef),digits=2),round(totdef,digits=2)]
   simresult[i,13:16] .= [round(sm,digits=4),round(ea,digits=4),round(outx,digits=4),round(outbog,digits=4)]
   simresult[i,17:20] .= [round(outglac,digits=2), round(m_r_onglac+outglac,digits=2),round(lyrs,digits=4), round(Qmm,digits=4)]
   simresult[i,21:24] .= [round(smbog,digits= 4), round(eabog,digits=4),round(boglyrs,digits=4),round(qmm_state,digits=4)]
   simresult[i,25:27] .= [round(wcs,digits=4), round(PotEvap,digits=4), round(MSWrad*fac, digits=4)]
   simresult[i,28:31] .= [round(MLT*fac,digits=4), round(MLA*fac,digits=4),round(MSH*fac,digits=4),round(MLE*fac,digits=4)]             
   simresult[i,32:35] .= [round(MALB,digits=4), round(midGLAC,digits=4), round(midGLAC1,digits = 4),round(mCl,digits = 4)]                
    else
   simresult[i, 1:5] .= [ptqinn.yr[i],ptqinn.mnt[i], ptqinn.day[i], ptqinn.hr[i],ptqinn.min[i] ]     
   simresult[i,6:9] .= [round(meanprecip,digits=4), round(meantemp,digits=3), ptqinn.q[i], round(Qm3s,digits=6)]
   simresult[i,10:13] .= [round(middelsca,digits=3),round(snomag,digits=4),round((M-totdef),digits=2),round(totdef,digits=2)]
   simresult[i,14:17] .= [round(sm,digits=4),round(ea,digits=4),round(outx,digits=4),round(outbog,digits=4)]
   simresult[i,18:21] .= [round(outglac,digits=2), round(m_r_onglac+outglac,digits=2),round(lyrs,digits=4), round(Qmm,digits=4)]
   simresult[i,22:25] .= [round(smbog,digits= 4), round(eabog,digits=4),round(boglyrs,digits=4),round(qmm_state,digits=4)]
   simresult[i,26:28] .= [round(wcs,digits=4), round(PotEvap,digits=4),round(MSWrad*fac, digits=4)] 
   simresult[i,29:32] .= [round(MLT*fac,digits=4), round(MLA*fac,digits=4),round(MSH*fac,digits=4),round(MLE*fac,digits=4)]             
   simresult[i,33:36] .= [round(MALB,digits=4),round(midGLAC,digits=4), round(midGLAC1,digits = 4),round(mCl,digits = 4)]               
    end
  end
  #grwpointresult[i, 1:29] <- c(ptqinn$yr[i],ptqinn$mnt[i],ptqinn$day[i],ptqinn$hr[i], round(GrWPoint[1:25],2))
  
    #end #Temporær end

  qberegn[i] = Qm3s      #Simulated runoff for current timestep. Used for skillscores etc.
          
  #---------------------------------------------------------------------------------------------------

#  #Saving states, SWE, sm, Layers, etc in a RData file
#  if(i == (savest[tell]) && savestate == 1) 
#  
#    Monthlast = ptqinn.mnt[savest[tell]]
#    Daylast = ptqinn.day[savest[tell]]
#    Yearlast = ptqinn.yr[savest[tell]]
#    Hourlast = ptqinn.hr[savest[tell]]
#   tdate <- paste(Yearlast,".",Monthlast,".",Daylast,".",Hourlast,sep="") 
#   tdate1 <- strptime(tdate, "%Y.%m.%d.%H")
#   tilsdate <- format(tdate1,"%Y%m%d.%H")
#   tilst.file2 <- paste(path2,"tilst\\","tilst_",i,"_",tilsdate,".rdata",sep="") 
#   tilst <-c("nsno","alfa","ny","sm","spd","wcd","sca", "Layers", "QRivx","smbog", "totdef","tilsdate","k")
#   save(list=tilst, file = tilst.file2)
#    if(tell < 1870)tell <- tell +1
#  end
  #---------------------------------------------------------------------------------------------
 
  
end        # end of loop for number of timesteps in timeseries

   println("nodaysLake=",nodaysLake)
 skillstart = startsim 
 days2 = days
 meansim = mean(qberegn[skillstart:days2])
 meanobs = mean(ptqinn.q[skillstart:days2])
    
# Computing skillscores NSE, kge, bias
 KGE = (kge(qberegn[skillstart:days2],ptqinn.q[skillstart:days2]))
 NSE = (nse(qberegn[skillstart:days2],ptqinn.q[skillstart:days2]))
 bias = (meansim/meanobs)

 #Recall, tprm is the parameter vector
                                    #pkr,      skr,      u,        pro,      TX,        CGLAC,     CritFlux,  rv
dfy = DataFrame(A=NSE,B=KGE,C=bias,D=tprm[1],E=tprm[2],F=tprm[3],G=tprm[4],H=tprm[5],II=tprm[6],JJ=tprm[7],KK=tprm[8],
    LL=tprm[9],MM=tprm[10], NN=tprm[11])
CSV.write(r2fil,DataFrame(dfy), delim = ';', header=false, append = true)

if (kal == 0)
    CSV.write(utfile,DataFrame(simresult,:auto), delim = ';', header=false)
end
  
#println(days)                
if(kal == 0)
 start2005 = 1
 #psum = sum(simresult[start2005:(start2005+days-1),5])#2005
 #println(psum)
 #easum = sum(simresult[start2005:(start2005+days-1),14])#2005
 #println(easum)
 #println("ratio sum_ea/sum_p =",round(easum/psum,digits= 3))
 println("M=",round(M,digits = 3))
end

return ptqinn.q, qberegn, KGE, NSE, bias

end # end of DDDModulFunc
