#DDD model Module based for hydrological predictions in urban areas
#
#----------------------------------------------------------------------------------
#     Script name:  DDD2018_Urban.R
#     Author:     Thomas Skaugen
#     Revised: 14.1.2020
#              
# Purpose   : DDD (Distance Distribution Dynamics) hydrological model simulates
#               1)saturated soil water flow, states of saturated (2-D)) and unsaturates zones
#               2) snow accumulation, distribution and melt
#               3) urban runoff for a given catchment.
#               4) saves the states for each day where precip > 0, i. e where SCHADEX might insert simulated extreme precip
#               5) Uses EB for snowmelt and evapotranspiration
#               6) estimates overland flow rivernetwork 
#               7) we will erform simulations for three reservoirs: 1)IP impermeable ares 2) P, permeable areas and 3) water in conduits (we need to know capacity)
#
# Model input: 1) catchment physical/geographical description, distance distributions for landscape features and river network
#              2) weather data, temperature,precipitation 
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

function DDDUrbanFunc(Gpar, startsim, tprm, prm ,ptqfile,utfile,r2fil, modstate,savestate, kal, spinuptime)


Lty = 2                      # number of landscape types in urban areas, for now Lty = 2, where Permeable (P) = 1 and ImPermeable (IP) =2, does not include wetlands

#-------------------------------------preprocess start------------------------------------------
#-----------------------------------------------------------------------------------------------

#####################VARIABLES FOR 10 ELEVATION ZONES and for Parts ###########################################
isoil = zeros(Lty,10)# precipitation and snowmelt from the elevation zones
swgt = zeros(Lty,10)# weights for input to soils pr each elevation zone NB bogs are not part of this yet
spd = zeros(Lty,10)
swe_h = zeros(Lty,10)   #SWE pr altitude level
wcd = zeros(Lty,10)
sca = zeros(Lty,10)
nsno = zeros(Lty,10)
alfa = zeros(Lty,10)
ny = zeros(Lty,10)
snowfree = zeros(Lty,10)  #Indicator variable 1 is sca  less than gca 0 if sca > gca
tempvec = zeros(Lty,10) # temperature vector  to save in statefile
PR = zeros(Lty,10)
PS = zeros(Lty,10)
MW = zeros(Lty,10)
MWGLAC = zeros(Lty,10)
SWrad = zeros(Lty,10)
LA = zeros(Lty,10)
LT = zeros(Lty,10)
Pa = zeros(10)

#Addition for Enegy balance snowmelt
snro = zeros(Lty,10)            #Snow density
melt = zeros(Lty,10)           #instead of degreeday melt (MW)
snowdepth = zeros(Lty,10)        #hmm..
sndens = zeros(Lty,10)           #hm ...
#################################Parameters as vectors#############################################
a0 = zeros(Lty) # scale parameter unit snow
n0 = zeros(Lty) # shape parameter unit snow
d = zeros(Lty)
Aprim = zeros(Lty)
taux = zeros(Lty)
MAD = zeros(Lty)          # really just one value, total for the entire catchment
CritFlux = zeros(Lty) # Critical volume/time unit to create overland flow channels
Ltyfrac = zeros(Lty+1) # areal fraction landscape type permeable area, includes forest and parks,fields, impermeable are built up areas and streets
Ltymax = zeros(Lty+1) #includes wetlands
Ltymid = zeros(Lty+1) #includes wetlands
Ltyz = zeros(Lty+1)   #includes wetlands
   
#################################States as vectors#################################################
M = zeros(Lty)
ICap = zeros(Lty)       # Infiltration capacity foor P[1] and IP[2] in [mm/s]
lyrs = zeros(Lty)
subsurface = zeros(Lty)
area = zeros(Lty+1)     # includes the area for wetlands
totdef = zeros(Lty)
Def = zeros(Lty)
PotEvap = zeros(Lty)
sm = zeros(Lty)
tosoil = zeros(Lty)
gmltosoil = zeros(Lty)
outx = zeros(Lty)
sm = zeros(Lty)
ea = zeros(Lty)
ea_sm = zeros(Lty)
ea_S = zeros(Lty)
OF = zeros(Lty)
Snowdepth = zeros(Lty)
Albedo = zeros(Lty)
snittT = zeros(Lty)
MLT = zeros(Lty)
MLA = zeros(Lty)
MSWrad = zeros(Lty)
snomag = zeros(Lty)
middelsca = zeros(Lty)
snofritt = zeros(Lty)
wcs = zeros(Lty)
scaobx = zeros(10)
hprecip = zeros(10)
htemp = zeros(10)
hfelt = zeros(10)

########################Reading parameters#########################################################
#Hyposgraphic curve and locations
hfelt[1] = prm.val[2]+(prm.val[3]-prm.val[2])/2.0
hfelt[2] = prm.val[3]+(prm.val[4]-prm.val[3])/2.0
hfelt[3] = prm.val[4]+(prm.val[5]-prm.val[4])/2.0
hfelt[4] = prm.val[5]+(prm.val[6]-prm.val[5])/2.0
hfelt[5] = prm.val[6]+(prm.val[7]-prm.val[6])/2.0
hfelt[6] = prm.val[7]+(prm.val[8]-prm.val[7])/2.0
hfelt[7] = prm.val[8]+(prm.val[9]-prm.val[8])/2.0
hfelt[8] = prm.val[9]+(prm.val[10]-prm.val[9])/2.0
hfelt[9] = prm.val[10]+(prm.val[11]-prm.val[10])/2.0
hfelt[10] = prm.val[11]+(prm.val[12]-prm.val[11])/2.0
hfeltmid = hfelt[5]             # mean elevation of catchment
phi = prm.val[13]*pi/180  #converts Lat degrees to radians lat
thi = prm.val[14]*pi/180  #converts Lon degrees to radians Lon 

#Meteorological corrections
pkorr = tprm[4] #prm.val[15]              # Met corrections rain
skorr = prm.val[16]               # Met corrections snow
u = tprm[1] #prm.val[17]                 # mean vind velocity [m/s]

#Snow parameters
pro = tprm[2] #prm.val[18]               # [fraction] liquid water content in snow
TX = tprm[3] #prm$V2[19]                 # Threshold temp for rain/snow
a0[1:Lty] .= prm.val[20:21]         # Alfa null in snofall unit distribution P, IP
d[1:Lty] .= prm.val[22:23]          # Decorrelation length (Measured in n) for precipitation P, IP
Aprim[1:Lty] .= prm.val[24:25]      # Initial albedo 0.86

taux[1:Lty] .= prm.val[26:27]       # Age of snow (UEB) 0.0
CX = 0.01  # degree-day factor for snowmelt Not operational
TS = 0.0  # thresholdtemperature for snowmelt Not operational

# Misc.
Timeresinsec = prm.val[28] # in seconds
MAD[1:Lty] .= prm.val[29:30]        # mean annual discharge
totarea = prm.val[31]            # Total area of catchment in m2
NoL = Int(prm.val[32])                # Number of layers in subsurface including OF level
R = prm.val[33]                  # soil moisture content/100 for field capacity of soils 
CritFlux[1:Lty] .= prm.val[34:35]  # Critical volume/time unit to create overland flow channels
#CritFlux[1] = tprm[4]
    
#Distance distributions Lty
Ltyfrac[1:(Lty+1)] .= prm.val[36:38]  # areal fraction permeable area:forest and parks,fields, impermeable area:roads and built up areas, wetlands 
Ltymax[1:(Lty+1)] .= prm.val[39:41]  # Maximum distances Landscapetypes
Ltymid[1:(Lty+1)] .= prm.val[42:44] 
#Ltymid[1] = tprm[5]  # mean distance permeable, Impermeable area, wetlands
Ltyz[1:(Lty+1)] .= prm.val[45:47]     # frac zero dist permable, impermeable area, wetlands

#Distance distribution conduit network
midFl = prm.val[48]        # mean of distance in river/conduits network
stdFl = prm.val[49]        # standard deviation of distances in river/conduits network
maxFl = prm.val[50]        # maximum distance in river/conduits network

#Celerities, subsurface and conduits
OFVP  = tprm[5]#prm.val[51]       # Overland flow velocity P
OFVIP = prm.val[52]       # Overland flow velocity IP
GshInt = tprm[6] #prm.val[53]        # shapeparameter saturation Capital lambda
GscInt = tprm[7] #prm.val[54]        # scaleparameter saturation Capital lambda
Gshape, Gscale = Big2SmallLambda(GshInt, GscInt) # Coverting integrated celerity to layers    

#Gshape = Gpar[1]
#Gscale = Gpar[2]          # Coverting integrated celerity to layers takes to long in calibration: preprocessing

rv = prm.val[55]    # celerity of water in river/conduits network, we'll probably want a celeritydistribution associated with levels here, have data!
persons = tprm[8] # prm.val[56]
    
#Hardcoded
#ICap[1] = (1000)*(GshInt*GscInt*Ltymid[1]/Timeresinsec) # Infiltration capacity for Permeable[1] [mm/s]. Equals mean Subsurface velocity
#ICap[2] = ICap[1]*0.001  # Infiltration capacity for ImPermeable[2] in

#provided by parameterfile
ICap[1] = prm.val[57]/3600.0#  Infiltration capacity for Permeable[1] Input from parameterfile is in mm/hour. we here transform to [mm/s]and later on to [mm/Timeresinsec] 
ICap[2] = prm.val[58]/3600.0# # Infiltration capacity for ImPermeable[2]  previously as ICap[1]*0.002  

#Infiltration capacity in mm/s. Must multiply with Timeresinsec
for Lst in 1:Lty 
    ICap[Lst] = ICap[Lst] .* Timeresinsec
end

#reading input file (ptq)
ptqinn = CSV.read(ptqfile, DataFrame, header=["yr","mnt","day","hr","min","p01","p02","p03","p04","p05","p06","p07","p08","p09","p10",
        "t01","t02","t03","t04","t05","t06","t07","t08","t09","t10","q"])    # reading ptq data, elevation zones
                                                                                                                                
days = Int(length(ptqinn.yr))    # Length of timeseries

#Constants
hson = 10                               # Elevation zones
unitsnow = 0.1                          # mean of unit snow [mm]
n0 = unitsnow*a0                        # shape parameter of unit snow
gtcel = 0.99                            # threshold groundwater capacity -> Overland flow
CFR = 2.5*(Timeresinsec/86400)*0.0833   # Fixed as 1/12 of estimate of CX= 2.5 for 24 hours 
len = Int(5*(86400/Timeresinsec))       # number of timestepes to spin up the model. Recommended to use timesteps equal to a minimum of 5 days to estimate the snowpack temperature.)
startsim = startsim + len               # taking into account estimating snowpack temperature
STempvec = zeros(len)              # Temperature vector for estimating snowpack temperature

Pa[1:10] = 101.3*((293 .- 0.0065.*hfelt[1:10])/293.0).^5.26     # Air pressure as a function of height Zhu et al. 2012 and Stoy et al, 2018
MPa = mean(Pa)

#vectors and matrixes
ddist = zeros(Lty,NoL)          # distribution of water in layers
ddistx = zeros(Lty,NoL)         # water present in layers
ddistxtemp = zeros(Lty,NoL)     # temporay of the above
Magkap = zeros(Lty,NoL)
aktMag = zeros(Lty,NoL)         # current water in Layers
k = zeros(Lty,NoL)              # subsurface celerities (including that of overland flow) 
qberegn = zeros(days)            # for calculation of skill scores etc. 

if(kal == 0)
    simresult = zeros(days,34)        # matrix which into the results are written
end 

#grwpointresult <-matrix(0.0, days,29)    # matrix which into the groundwater point results are written

# Areas and weights
area = Ltyfrac.* totarea                   #areas of differnt Lty, a vector

elevarea = (totarea/hson)               #area pr elevationzone
    

#Finds the fraction of soils (and bogs) in each elevation zone in relation to total soil (and bog) area

for Lst in 1:Lty
   swgt[Lst,1:hson] .= elevarea/totarea           #weigths relative to total area of landscape types pr elevation area
end
                  
#Subsurface celerities 

for Lst in 1:Lty
    k[Lst,1:NoL] = CeleritySubSurface(NoL, Gshape, Gscale, Ltymid[Lst], Timeresinsec) # Celerity of subsurface (and overland) flow
end

#Overlandflow celerities
k[1,1] = OFVP  # P m/s Holden et al. WRR,2008,  # Overland flow celerites P
k[2,1] = OFVIP # IP, Sedyowati et al. 2017      # Overland flow celerites IP#Overlandflow celerities

bogspeed = k[1,1] # same as for overland flow Permeable areas

                                            
#Conduit celerities/velocities

#Unit hydrgraphs for Wetlands
antBogsteps =Int(trunc(Ltymax[3]/(bogspeed*Timeresinsec))+1) # Number of timsteps to drain wetlands
UHbog = zeros(antBogsteps)                  # Vector for Bog unit hydrograph
if(area[3] > 0)
    UHbog = SingleUH(bogspeed,Timeresinsec, Ltymid[3], Ltymax[3], Ltyz[3])
end
if(area[3] <= 0)
    UHbog = 1.0
end

#Unit hydrographs for P based on celerities and distance distributions,exponentially distributed
antHorlag = zeros(Int64,Lty)
for Lst in 1:Lty
    antHorlag[Lst] = Int(trunc(Ltymax[Lst]/(k[Lst,NoL]*Timeresinsec))+1)       # +1 includes the final day.  NB only in soils part
end

nodaysvector = zeros(Int64,Lty, NoL)                        # Number of timesteps to drain each layer
for Lst in 1:Lty 
    nodaysvector[Lst,1:NoL] = Int.(trunc.(Ltymax[Lst]./k[Lst,1:NoL]./Timeresinsec).+1)   # integer time steps
end


layerUH_P = zeros(NoL,antHorlag[1])
layerUH_IP = zeros(NoL,antHorlag[2])  

for i in 1: NoL  
  for Lst in 1:1  
   layerUH_P[i,1:nodaysvector[Lst,i]] .= SingleUH(k[Lst,i], Timeresinsec, Ltymid[Lst], Ltymax[Lst], Ltyz[Lst])
  end
  for Lst in 2:2 
   layerUH_IP[i,1:nodaysvector[Lst,i]] .= SingleUH(k[Lst,i], Timeresinsec, Ltymid[Lst], Ltymax[Lst], Ltyz[Lst])
  end
end

#UH for Conduits, normal distributed 
#number of time units in river routing, rounded up. If set equal to one day, time is actually less
nodaysRiv = Int(trunc(maxFl/rv/Timeresinsec)+1)
println(nodaysRiv)
if (nodaysRiv > 1)
  timeres = [0 : 1: (nodaysRiv-1);]                      # For UH, same as antbox.
  midFlscl = midFl/rv/Timeresinsec
  stdFlscl = stdFl/rv/Timeresinsec
  UHriver = zeros(nodaysRiv)
  Rivnorm = Normal(midFlscl,stdFlscl)
  UHriver = pdf.(Rivnorm,timeres)                       # makes  the PDF
  UHriver .= UHriver/sum(UHriver)                         # scale to give unity sum = 1
  noDT = Int(length(UHriver))
else
  UHriver = 1.0                      # implies no river routing
  noDT = 1
end    

qsimutx = zeros(noDT)      # common vector for runoff from all landscapetypes [m3/s]
qsimutxP = zeros(noDT)     # vector for runoff permeabel i m3/s
qsimutxIP = zeros(noDT)    # vector for runoff impermeable i m3/s
qsimutxBog = zeros(noDT)   # vector for runoff from wetlands bogss
qsimutxOF = zeros(noDT)    # vector for Overland Flow
QRivx = zeros(noDT)        # all qsimutvectors are stored UHriver timesteps back. 
QRivxP = zeros(noDT)
QRivxIP = zeros(noDT)
QRivxBog = zeros(noDT)
QRivxOF = zeros(noDT)

#Groundwater layers; 2dim levels, 1 fastest, NoL slowest#
LayersP = zeros(NoL,antHorlag[1])
LayersIP = zeros(NoL,antHorlag[2])

BogLayers = zeros(antBogsteps)

for Lst in 1: Lty
  Magkap[Lst,1:NoL], M[Lst] = LayerEstimation(GshInt,GscInt,Timeresinsec,Ltymax[Lst],Ltymid[Lst], 
        Ltyfrac[Lst]*MAD[Lst], Ltyfrac[Lst]*totarea, NoL, gtcel) 
end
#Note tha MAD is for the entire catchment. It is fractioned in the above subroutine

#Extracting areas and heights in pyramide-plot
AreasP, delta_dP = PyrAreas(NoL,Ltyfrac[1]*totarea,Ltymax[1],nodaysvector[1,1:NoL], layerUH_P, antHorlag[1]) # A matrix of areas[in m2] for each time-step box, 
                                                                 # Important: same dimensions as Layers
                                                                 # The area drained pr time-step box for each layer

AreasIP, delta_dIP = PyrAreas(NoL,Ltyfrac[2]*totarea,Ltymax[2],nodaysvector[2,1:NoL], layerUH_IP, antHorlag[2]) # A matrix of areas[in m2] for each time-step box, 

#States
modstate = 0
#if(modstate ==1) load(tilst.file)  # Run with state
if (modstate == 0)                   # Run with initial values  
   smbog  = 0.0                          # Soilmoisture Bogs   
   QRivx[1:noDT] .= qsimutxP[1:noDT]             #discharge in m3/s, contribution from P
   QRivx[1:noDT] .= QRivx[1:noDT] .+ qsimutxIP[1:noDT] #adding contribution from IP
   QRivx[1:noDT] .= QRivx[1:noDT] .+ qsimutxBog[1:noDT]#adding contribution from bogs
   QRivxBog[1:noDT] .= qsimutxBog[1:noDT] 
   QRivxP[1:noDT] .= qsimutxP[1:noDT] 
   QRivxIP[1:noDT] .= qsimutxIP[1:noDT]
   totdef[1:Lty] .= M[1:Lty]
end

Snowdepth[1:Lty] .= 0.0
snittT[1:Lty] .= 0.0 


############################################################################
#                                                                          #
#simulation starts here.............. days is no. time steps               #
############################################################################
for i in startsim:days


  dato = Date(ptqinn.yr[i],ptqinn.mnt[i], ptqinn.day[i]) #date
  #println(i," ",dato)      
        
  hr = ptqinn.hr[i] 
  DN = Dates.dayofyear(dato)           #daynumber        
  
  #Reads Precipitation and temperature for each elevation zone   
  htemp = ptqinn[i,16:25]
  hprecip = ptqinn[i,6:15] 

  meanprecip = mean(hprecip)
  meantemp =  mean(htemp)
  
  for Lst in 1:Lty     # landscape types, one snow regime for each landscape type

    for idim in 1:hson # elevation zones  

      STempvec = TemperatureVector(ptqinn[(i-len+1):i,16:25],idim, len)
      CGLAC = 0.0 # dummy, has no role

      #Calculates Rain vs Snow and Snow/Glaciermelt
      WaterIn = NedbEB(CFR,DN,hprecip[idim],htemp[idim],spd[Lst,idim],Timeresinsec,TX,pkorr,skorr,
                         hr,u,Pa[idim],CGLAC,CX,TS,taux[Lst],snittT[Lst],phi,thi)

      #From WaterIn per elevation zone:         
      PS[Lst,idim] = WaterIn[1]         # Precipitation as snow
      PR[Lst,idim] = WaterIn[2]         # Precipitation as rain      
      MW[Lst,idim] = WaterIn[3]         # Snow melt potential
      #MWGlac[LST,idim] = WaterIn[4]  # glaciers, not relevant
      SWrad[Lst,idim] = WaterIn[5]    # Short wave radiation (net)
      LA[Lst,idim] = WaterIn[6]         # Long wave radiation atomspheric
      LT[Lst,idim] = WaterIn[7]          # Long wave radiation terrestrial
      SH = WaterIn[8]          # Turbulent heat
      LE = WaterIn[9]          # Latent heat
      GH = WaterIn[10]          # Ground heat
      PH = WaterIn[11]          # Precipitation heat
      CC = WaterIn[12]          # Cold content
      Albedo[Lst] = WaterIn[13]            # Albedo
      taux[Lst] = WaterIn[14]     # Aging factor of snow for albedo calulations
      snittT[Lst] = WaterIn[15]  # Snowpack temperature
      Tss = WaterIn[16]        # Snow surface temperature
      RH =  WaterIn[18]         # Relative humidity
      Cl = WaterIn[17]          # CloudCover

      # calculates Snow accumulation and Snowdistribution     
      FromSnow = SnowGamma(PR[Lst,idim],PS[Lst,idim],MW[Lst,idim],sca[Lst,idim],spd[Lst,idim],wcd[Lst,idim],
                        pro,nsno[Lst,idim],alfa[Lst,idim],ny[Lst,idim],a0[Lst],n0[Lst],a0[Lst],d[Lst])
    
      isoil[Lst,idim] = FromSnow[1]    #precipitation and snowmelt
      spd[Lst,idim] = FromSnow[2]    # snow water equivalent
      wcd[Lst,idim] = FromSnow[3]    # liquid water in snow
      sca[Lst,idim] = FromSnow[4] 
      nsno[Lst,idim] = FromSnow[5] 
      alfa[Lst,idim] = FromSnow[6] 
      ny[Lst,idim] = FromSnow[7]
     
      if (spd[Lst,idim] > 20000)     
        spd[Lst,idim] = 20000.0
        println("Glaciers are made. This speeds up calibration")
      end
        
     #Snowdepths and snowdenities
     snowdepthx = snowdepth[Lst,idim]
     snomag[Lst] =  sca[Lst, idim]*spd[Lst,idim] #snowreservoir this timestep

     #Calling snowdepth routine
     snowdepth[Lst,idim],sndens[Lst,idim] = NewSnowSDEB(PR[Lst,idim],PS[Lst,idim],htemp[idim],TX,
            MW[Lst,idim],snomag[Lst],snomag[Lst],snowdepthx)
     
    end #for elevation zones
 
    MSWrad[Lst] = mean(SWrad[Lst,1:hson])
    MLA[Lst] = mean(LA[Lst,1:hson])
    MLT[Lst] = mean(LT[Lst,1:hson])

    #Snowreservoir
    snomag[Lst] = mean(sca[Lst,1:hson] .* spd[Lst,1:hson]) #mean catchment SWE ,must multiply with SCA to have arealvalues
    swe_h[Lst,1:hson] = sca[Lst,1:hson] .* spd[Lst,1:hson] #SWE pr. elevation zone
    middelsca[Lst] = sum(sca[Lst,1:hson])/hson
    snofritt[Lst] = 1-middelsca[Lst]
    wcs[Lst] = mean(wcd[Lst,1:hson])
  

  end #For landscape types snow and rain
    

    #Estimating current capacity in Layers before Evapotranspiration
    ddistx[1,1:NoL] = LayerCapacityUpdate(LayersP, nodaysvector[1,1:NoL], Magkap[1,1:NoL], NoL)
    ddistx[2,1:NoL] = LayerCapacityUpdate(LayersIP, nodaysvector[2,1:NoL], Magkap[2,1:NoL], NoL)
   
    #Updating the deficit (for all sub surface layers, NOT overland flow layer) before Evapotranspiration
    for Lst in 1:Lty
      totdef[Lst] = sum(ddistx[Lst,2:NoL])            # This may be negative if input outstrips deficits
      gmltosoil[Lst] = sum(isoil[Lst,1:hson] .* swgt[Lst,1:hson])  #snowmelt and rain on soils. isoil is weighted by the fraction of soils. Is later multiplied with area for the different fractions.
    
      PotEvap[Lst] = PotentialEvapPT(meantemp, MSWrad[Lst], MLA[Lst], MLT[Lst], MPa)
      if(PotEvap[Lst] < 0.0)
        PotEvap[Lst] = 0.0
      end
      Def[Lst] = max(totdef[Lst],0)
    end
	
	# Assigning overland flow if infiltration capapcity is exceeded and reduces tosoil in the process
	for Lst in 1:Lty
      OF[Lst],tosoil[Lst] = OFICap(gmltosoil[Lst],ICap[Lst])
	  #if (meanprecip > 0.0)
	   # println("precip=", meanprecip)
       # println("OF=",OF)
	   # println("tosoil=", tosoil)
	   # println(i," ",dato)
	  #end
    end 	  
   
    #calculating evapotranspiration from soilmoisture reservoir
    for Lst in 1:Lty
     sm[Lst], ea_sm[Lst], ea_S[Lst] = UnsaturatedEvapEB(tosoil[Lst], meantemp, sm[Lst], M[Lst], Def[Lst], PotEvap[Lst])
                  # updated sm
                  # evapotranspiration to be drawn from the soilmoisture
                  # evapotranspiration to be drawn from Layers  
    end
    
    #calculating additional (ea_S) evapotranspiration directly from Layers
    LayersP, ea_S[1] = LayerEvap(LayersP, nodaysvector[1,1:NoL], ea_S[1],layerUH_P,NoL)
    LayersIP, ea_S[2] = LayerEvap(LayersIP, nodaysvector[2,1:NoL], ea_S[2],layerUH_IP, NoL)
    ea[1] = ea_sm[1] + ea_S[1]
    ea[2] = ea_sm[2] + ea_S[2] # adding evapotranspiration from sm (first) and Layers (second)
       
    #Estimating current capacity in Layers after Evapotranspiration
    ddistx[1,1:NoL] = LayerCapacityUpdate(LayersP, nodaysvector[1,1:NoL], Magkap[1,1:NoL], NoL)
    ddistx[2,1:NoL] = LayerCapacityUpdate(LayersIP, nodaysvector[2,1:NoL], Magkap[2,1:NoL], NoL)
        
    #Updating the deficit (for all sub surface layers, NOT overland flow layer)
    for Lst in 1:Lty       
      totdef[Lst] = sum(ddistx[Lst,2:NoL])            # This may be negative if input outstrips deficits
      Def[Lst] = max(totdef[Lst],0)
    end 
     
  #Estimate isoil from all elevation zones and Lty
  for Lst in 1:Lty   
   # State of unsaturated zone, it can etiher be zero (complete saturation) or positive, the deficit. If negative we have overland flow 
   smlast = sm[Lst]
   # Call for the soilwater routine  calculates soil water (unsaturated) and evapotransipration 
   outx[Lst], sm[Lst] = UnsaturatedExEvap(tosoil[Lst], sm[Lst], R, Def[Lst])

   #if(round((tosoil[Lst]),digits=3)-round((outx[Lst]+(sm[Lst]-smlast)),digits = 3)!= 0.0)  # balance OK 20.06.2023     
   #  println(i,"hi, tosoil and sm not in balance", Lst)
   #  println(tosoil[Lst])
   # println(outx[Lst]+(sm[Lst]-smlast))
   #end   
  end
    
  #From Wetlands/Bogs
   smboglast = smbog
   Fromwetlands = WetlandsEB(tosoil[1], meantemp, middelsca[1], smbog, M[1], PotEvap[1])# takes on Permeable values, Lst =1
   outbog = Fromwetlands[1]
   smbog = Fromwetlands[2]
   eabog = Fromwetlands[3]

   #if(round((tosoil[1]),digits= 3)-round((outbog+(smbog+eabog-smboglast)),digits=3)!= 0.0)  # balance so far   
   #  println(i," hei after WetlandsEB")
   #  println(tosoil[1])
   #  println(outbog+(smbog+eabog-smboglast))
   #  println(sum(BogLayers))  #Boglayers is the routing of runoff generated from bogs, NOT the reservoir
   # end
  
  #Estimating input to Layers according to capacity
  for Lst in 1:Lty  
   ddist[Lst,1:NoL] = GrvInputDistributionICap(outx[Lst], NoL, ddistx[Lst,1:NoL], ddist[Lst,1:NoL], ICap[Lst])
  end
  
  #We need to include OF due to exceeeded infiltration capacity and UPDATE dist and outx. ddist will distribute outx taking into account both sat and inf excess
  for Lst in 1:Lty  
   if(OF[Lst]>0.0) # recall OF is only due to inf excess
        ddist[Lst,1] = ddist[Lst,1] + (OF[Lst]/(OF[Lst]+outx[Lst])) # will include water from exceedance of infiltration capacity
		ddist[Lst,2:NoL] = (1-ddist[Lst,1]).*ddist[Lst,2:NoL]  #updating the rest of ddist if OF
	    outx[Lst]= outx[Lst]+ OF[Lst]
   end
  end
  
  #Calculating UH for overland flow Layer #1
  OFDDD = 1
  if (OFDDD ==1)  
    for Lst in 1:Lty
     if(ddist[Lst,1]*outx[Lst] > 0)# overland flow for LST =1 (P) is confirmed for this event
      
       if(Lst==1)
        layerUH_P[1, 1:nodaysvector[1,1]] = OverlandFlowDynamicDD(k[Lst,1:NoL],ddist[Lst,1:NoL],outx[Lst], layerUH_P,
                    nodaysvector[Lst,1:NoL],NoL, Ltymid[Lst], CritFlux[Lst], Timeresinsec)
       end              
       if(Lst==2)
        layerUH_IP[1, 1:nodaysvector[2,1]] = OverlandFlowDynamicDD(k[Lst,],ddist[Lst,],outx[Lst], layerUH_IP,
                    nodaysvector[Lst,1:NoL],NoL, Ltymid[Lst], CritFlux[Lst], Timeresinsec)
       end                
       
     end
    end
  end

  ##Waterbalance calculations
  #RinnP =sum(QRivxP[2:noDT])*(Timeresinsec*1000/area[1]) # P water stored in RN from last TimeSep
  #RinnIP= sum(QRivxIP[2:noDT])*(Timeresinsec*1000/area[2]) # IP water stored in RN from last TimeStep
        
  #SPinn = sum(LayersP) + sum(ddist[1,1:NoL] .* outx[1] .* layerUH_P[1:NoL,1:antHorlag[1]])#Alle lag i P pluss dagens input før runoff
  ##SPinn1 = sum(LayersP[1,1:nodaysvector[1,1]]) + sum(ddist[1,1] .* outx[1] .* layerUH_P[1,1:nodaysvector[1,1]])#lag1
  ##SPinn2 = sum(LayersP[2,1:nodaysvector[1,2]]) + sum(ddist[1,2] .* outx[1] .* layerUH_P[2,1:nodaysvector[1,2]])#lag2 
  #SIPinn = sum(LayersIP) + sum(ddist[2,1:NoL] .* outx[2] .* layerUH_IP[1:NoL,1:antHorlag[2]])#Alle lag i IP pluss dagens input før runoff
  ##SIPinn1 = sum(LayersIP[1,1:nodaysvector[2,1]]) + sum(ddist[2,1] .* outx[2] .* layerUH_IP[1,1:nodaysvector[2,1]])
  ##SIPinn2  = sum(LayersIP[2,1:nodaysvector[2,2]]) + sum(ddist[2,2] .* outx[2] .* layerUH_IP[2,1:nodaysvector[2,2]])

    
  GDT_P = sum(LayersP[1:NoL,1]) + sum(ddist[1,1:NoL] .* outx[1] .* layerUH_P[1:NoL,1])      #groundwater to be discharged into the rivernetwork + this timesteps contribution
  GDT_IP = sum(LayersIP[1:NoL,1]) + sum(ddist[2,1:NoL] .* outx[2] .* layerUH_IP[1:NoL,1])   #groundwater to be discharged into the rivernetwork 
  GDT_Bog =  BogLayers[1]+outbog*UHbog[1]         #bogwater to be discharged into the rivernetwork + this timesteps contribution
  GDT_OF = Ltyfrac[1]*(LayersP[1,1]+(ddist[1,1]*outx[1]*layerUH_P[1,1]))+
           Ltyfrac[2]*(LayersIP[1,1]+(ddist[2,1]*outx[2]*layerUH_IP[1,1]))+
           Ltyfrac[3]*(BogLayers[1]+outbog*UHbog[1])# Overland flow part 

  #Updating the saturation Layers
    for Lst in 1:Lty
      if(Lst==1)
         LayerUpdate!(ddist[Lst,1:NoL],outx[Lst], LayersP, layerUH_P, nodaysvector[Lst,1:NoL], NoL)
    end       
      if(Lst==2)
         LayerUpdate!(ddist[Lst,1:NoL],outx[Lst], LayersIP, layerUH_IP, nodaysvector[Lst,1:NoL], NoL)       
      end
    end  

  BogLayerUpdate!(outbog, BogLayers, UHbog, antBogsteps)#
   
  #summing up groundwater states
    for Lst in 1:Lty
      if(Lst == 1)
         lyrs[Lst] = sum(LayersP)
         subsurface[Lst] = sum(LayersP[2:5,1:antHorlag[1]]) # gives the sum of layers after todays runoff has taken place
      end
      if(Lst==2)
         lyrs[Lst] = sum(LayersIP)
         subsurface[Lst] = sum(LayersIP[2:5,1:antHorlag[2]])
      end
    end
  
  #summing up bogstates
  boglyrs = sum(BogLayers)
        
  #reinstate original overland flow layer          
  for Lst in 1:1  
   layerUH_P[1,1:nodaysvector[Lst,1]] .= SingleUH(k[Lst,1], Timeresinsec, Ltymid[Lst], Ltymax[Lst], Ltyz[Lst])
  end
  for Lst in 2:2 
   layerUH_IP[1,1:nodaysvector[Lst,1]] .= SingleUH(k[Lst,1], Timeresinsec, Ltymid[Lst], Ltymax[Lst], Ltyz[Lst])
  end
  
  #runoff
  qsimutxP .= (((GDT_P/1000)*area[1])/Timeresinsec) .* UHriver     #[m3/s]
  qsimutxIP .= (((GDT_IP/1000)*area[2])/Timeresinsec) .* UHriver   #[m3/s]
  qsimutxBog .= (((GDT_Bog/1000)*area[3])/Timeresinsec) .* UHriver #[m3/s]
  qsimutxOF .= (((GDT_OF/1000)*totarea)/Timeresinsec) .* UHriver   #[m3/s]  # Overland flow. Not a contribution, just extracted the fraction OF of total runoff
  
  #Total response
  qsimutx[1:noDT] .= 0.0                                   # this is written anew for each timestep
  qsimutx[1:noDT] .= qsimutx[1:noDT] .+ qsimutxP[1:noDT]   #adding contribution from permeable areas
  qsimutx[1:noDT] .= qsimutx[1:noDT] .+ qsimutxIP[1:noDT]  #adding contribution from impermeable areas
  qsimutx[1:noDT] .= qsimutx[1:noDT] .+ qsimutxBog[1:noDT] #adding contribution from bogs

  
  #Updating the routing in river network total)
  QRivx, QRD = RiverUpdate(noDT,QRivx,qsimutx)
  qmm_state = (sum(QRivx)-QRD)*(Timeresinsec*1000/totarea)       # This is also a reservoir [mm], water from todays event is stored for future runoff in the RN, relevant for WB.
                                                                  # Minus QRD since this will be handled in the WB equations
  #Updating the routing in river network Bogs
  QRivxBog,QRDBog = RiverUpdate(noDT,QRivxBog,qsimutxBog)
    
  #Updating the routing in river network P
  QRivxP,QRDP = RiverUpdate(noDT,QRivxP,qsimutxP)
    
  #Updating the routing in river network IP
  QRivxIP,QRDIP = RiverUpdate(noDT,QRivxIP,qsimutxIP)
  
  #Updating the routing in river network OF
  QRivxOF,QRDOF = RiverUpdate(noDT,QRivxOF,qsimutxOF)
  
  Qmm = (QRD*Timeresinsec*1000/totarea)  #QRD in mm/day

   #WBP = (SPinn + RinnP) - (lyrs[1] + sum(QRivxP)*(Timeresinsec*1000/area[1]))    #mm  last term inludes discharge and storage in rivernetwork                                         
   #println("WBP = ", WBP, " Total inn P= ",SPinn)                                             
   #WBIP = (SIPinn + RinnIP) - (lyrs[2] + sum(QRivxIP)*(Timeresinsec*1000/area[2]))   #mm  last term inludes discharge and storage in rivernetwork                                         
   #println("WBIP = ", WBIP," Total inn IP= ",SIPinn)  
##  Water balance before and after runoff calculations ok 20.06.2023. Also check that precip-(evap + sm + outx) == 0

  if (kal==0)
  #Assigning outdata to vector, one vector for each timestep
   simresult[i, 1:5] = [ptqinn.yr[i],ptqinn.mnt[i],ptqinn.day[i],ptqinn.hr[i],ptqinn.min[i]]
   simresult[i, 6:9] = [round(meanprecip,digits= 3),round(meantemp,digits = 3), ptqinn.q[i],round(QRD,digits=6)]
   simresult[i, 10:13] = [round(QRDP,digits = 6),round(QRDIP,digits = 6), round(QRDBog,digits = 6),round(middelsca[1],digits=3)]
   simresult[i, 14:16] = [round(snomag[1],digits = 6), round(lyrs[1],digits = 6),round(lyrs[2],digits = 6)]
   simresult[i, 17:20] = [round(totdef[1],digits=2), round(totdef[2],digits=2),round(sm[1],digits=6),round(sm[2],digits=6)]
   simresult[i, 21:24] = [round(ea[1],digits=6),round(ea[2],digits=6),round(Qmm,digits=6),round(smbog,digits=6)]
   simresult[i, 25:28] = [round(eabog,digits=6),round(qmm_state,digits=6), round(boglyrs,digits=6), round(middelsca[2],digits=3)]
   simresult[i, 29:32] = [round(snomag[2],digits= 6), round(wcs[1],digits=6),round(wcs[2],digits=6),round(subsurface[1],digits=6)]
   simresult[i, 33:34] = [round(subsurface[2],digits=6),round(QRDOF,digits=6)]
  end
  #grwpointresult[i, 1:29] <- c(ptqinn$yr[i],ptqinn$mnt[i],ptqinn$day[i],ptqinn$hr[i], round(GrWPoint[1:25],2))
                                
  qberegn[i] = QRD      #Simulated runoff for current timestep. Used for skillscores etc.

  #---------------------------------------------------------------------------------------------------

  #Saving states, SWE, sm, Layers, etc in a RData file
  
 # if(i == (savest[tell]) && savestate == 1) 
 # {
 #   Monthlast<-ptqinn$mnt[savest[tell]]
 #   Daylast<-ptqinn$day[savest[tell]]
 #   Yearlast<-ptqinn$yr[savest[tell]]
 #   Hourlast<-ptqinn$hr[savest[tell]]
 #   tdate <- paste(Yearlast,".",Monthlast,".",Daylast,".",Hourlast,sep="") 
 #   tdate1 <- strptime(tdate, "%Y.%m.%d.%H")
 #   tilsdate <- format(tdate1,"%Y%m%d.%H")
 #   tilst.file2 <- paste(path2,"tilst\\","tilst_",i,"_",tilsdate,".rdata",sep="") 
 #   tilst <-c("nsno","alfa","ny","sm","spd","wcd","sca", "Layers", "QRivx","smbog", "totdef","tilsdate","k")
 #   save(list=tilst, file = tilst.file2)
 #    if(tell < 1870)tell <- tell +1
 # }
  #---------------------------------------------------------------------------------------------

end # end of loop for number of timesteps in timeseries
 
 skillstart = Int(spinuptime*(86400/Timeresinsec))
 days2 = days
 
wwater = 1.0*persons*(140.0 + 42.0)/(86400.0*1000.0)# norsk vann equivalent use in L/day add 30% leakage from water net

 meansim = mean(wwater .+ qberegn[skillstart:days2])
 meanobs = mean(ptqinn.q[skillstart:days2])

# Computing skillscores NSE, kge, bias
 KGE = (kge((wwater .+ qberegn[skillstart:days2]),ptqinn.q[skillstart:days2]))
 NSE = (nse((wwater .+ qberegn[skillstart:days2]),ptqinn.q[skillstart:days2]))
 bias = (meansim/meanobs)

#Recall, tprm is the parameter vector
#                        u,        pro,         TX,      PCritFlux,         OFP,        GshInt,   GscInt, persons
dfy = DataFrame(A=NSE,B=KGE,C=bias,D=tprm[1],E=tprm[2],F=tprm[3],G=tprm[4], H=tprm[5], II =tprm[6],IJ=tprm[7],
    JJ=tprm[8])
CSV.write(r2fil,DataFrame(dfy), delim = ';', header=false, append = true)

 if (kal==0)
  CSV.write(utfile,DataFrame(simresult,:auto), delim = ';', header=false)                
  println("M= ", M[1],M[2])
  println("k_P= ", k[1,1], k[1,2],k[1,3],k[1,4],k[1,5])
  println("k_IP= ", k[2,1], k[2,2],k[2,3],k[2,4],k[2,5])
  println("Mean(Qsim)= ",meansim)
 end           

return ptqinn.q, qberegn, KGE,NSE,bias

end  #end of func
