#Runs DDD model either in single run model or for calibration
#The model itself is called as a function which calls on several functions
#This one is for ECCO use not MRT and not operative or HSO
using CSV
using Distributions
using LsqFit
using Statistics
using Dates
using DataFrames
using Plots
using BlackBoxOptim
using JLD2

# Preprocessing routines
path_dddfun = joinpath(dirname(@__DIR__), "DDDFunctions")
include(joinpath(path_dddfun, "Big2SmallLambda.jl"))
#include(joinpath(path_dddfun, "CeleritySubSurface_DupBous.jl"))
include(joinpath(path_dddfun, "CeleritySubSurface.jl"))
include(joinpath(path_dddfun, "SingleUH.jl"))
include(joinpath(path_dddfun, "SingleNormalUH.jl"))
include(joinpath(path_dddfun, "LayerEstimation.jl"))
include(joinpath(path_dddfun, "PyrAreas.jl"))
include(joinpath(path_dddfun, "GrWPoint.jl"))
include(joinpath(path_dddfun, "RiverPoint.jl"))
include(joinpath(path_dddfun, "TemperatureVector.jl"))

##EB and Snow Routines
include(joinpath(path_dddfun, "NedbEBGlac_debug04072022.jl"))
include(joinpath(path_dddfun, "SnowpackTemp.jl"))
include(joinpath(path_dddfun, "TempstartUpdate.jl"))
include(joinpath(path_dddfun, "SmeltEBGlac_debug04072022.jl"))
include(joinpath(path_dddfun, "CloudCoverGlac_debug04072022.jl"))
include(joinpath(path_dddfun, "TssDewpoint.jl"))
include(joinpath(path_dddfun, "SolradTransAlbedoper_hrs_debug04072022.jl"))
include(joinpath(path_dddfun, "LongWaveRad_debug04072022.jl"))
include(joinpath(path_dddfun, "SensibleLatHeat_debug04072022.jl"))
include(joinpath(path_dddfun, "AlbedoUEB_debug04072022.jl"))
include(joinpath(path_dddfun, "GroundPrecCC.jl"))
include(joinpath(path_dddfun, "SnowGamma.jl"))
include(joinpath(path_dddfun, "Varc.jl"))
include(joinpath(path_dddfun, "NewSnowDensityEB.jl"))
include(joinpath(path_dddfun, "NewSnowSDEB.jl"))
include(joinpath(path_dddfun, "DensityAge.jl"))

#Subsurface and Evaporation routines
include(joinpath(path_dddfun, "LayerCapacityUpdate.jl"))
include(joinpath(path_dddfun, "PotentialEvapPT.jl"))
include(joinpath(path_dddfun, "UnsaturatedEvapEB.jl"))
include(joinpath(path_dddfun, "LayerEvap.jl"))
include(joinpath(path_dddfun, "UnsaturatedExEvap.jl"))
include(joinpath(path_dddfun, "WetlandsEB.jl"))
include(joinpath(path_dddfun, "GrvInputDistributionICap2022.jl"))
include(joinpath(path_dddfun, "OFICap.jl"))
include(joinpath(path_dddfun, "LayerUpdate.jl"))
include(joinpath(path_dddfun, "BogLayerUpdate.jl"))
include(joinpath(path_dddfun, "RiverUpdate.jl"))
## Overland Flow routine
include(joinpath(path_dddfun, "OverlandFlowDynamicDD.jl"))
## Efficiency criteria
include(joinpath(path_dddfun, "NSE_ths.jl"))
include(joinpath(path_dddfun, "KGE_ths.jl"))
# Model Module
#include(joinpath(path_dddfun, "DDDUrbanFunc.jl"))
include(joinpath(path_dddfun, "DDDAllTerrain22012024.jl")) #This one for ECCO use not MRT
########################################################################################

# Functions
function calib_wrapper_model(Gpar,startsim, tprm, prm, ptqfile, utfile, r2fil, modstate, savestate, kal, spinup)
 qobs, qberegn, KGE, NSE, bias = DDDAllTerrain(Gpar,startsim, tprm, prm, ptqfile, utfile, r2fil, modstate, savestate,
        kal, spinup)  
 return qobs,qberegn, KGE,NSE,bias 
end

function calib_single_wsh(Gpar,startsim, tprm, prm, ptqfile, utfile, r2fil, modstate, savestate, kal, spinup)
 qobs, qberegn, KGE, NSE, bias = DDDAllTerrain(Gpar,startsim, tprm, prm, ptqfile, utfile, r2fil, modstate, savestate,
        kal, spinup)    
 return (1.0 - KGE)
end

# Temporary data folders
dir_input = "/hdata/grid/DDD_test/inndataV2"
dir_param = "/hdata/grid/DDD_test/DDDurbanParameters"
root_output = "/hdata/grid/DDD_test/output"

# Options
catchment = "55.4"  # stationnumber
dir_output = mkpath(joinpath(root_output, catchment))
startsim = 1 
kal = 1
modstate = 0
savestate = 0
spinup = (31*4) #days used to spin up the model. 
resolution_time = "3h"         # this is just a marker for naming files, does NOT set the temporal resolution
ptqfile = joinpath(dir_input, catchment, string(catchment, "_" , resolution_time, "_ptq_DDDv2_kal.csv"))
paramfile = joinpath(dir_param, catchment, string("ParDDDv2_", catchment, "_3h.csv"))
r2fil = joinpath(dir_output, string("r2_", resolution_time, ".csv"))
utfile = joinpath(dir_output, string("simres_", resolution_time, ".csv"))

# Read parameter file
prm = CSV.read(paramfile, DataFrame, header=["Name", "val"], delim=';')
#            u,          pro           TX,         Pkorr        skorr,     GscInt      OVP          
#         OVIP       Lv            rv        
tprm = [prm.val[20], prm.val[21], prm.val[22], prm.val[18], prm.val[19],prm.val[33], prm.val[34], 
        prm.val[35],prm.val[36],prm.val[37]]
println(tprm)
Gshape, Gscale = Big2SmallLambda(prm.val[32], prm.val[33]) # Coverting integrated celerity to layers takes too long in calibration: preprocessing
Gpar = [Gshape, Gscale]
println(prm.val[32]," ", prm.val[33])


# Run or calibrate model
t1 = time()
if(kal == 0)
    qobs,qberegn,KGE,NSE, bias = calib_wrapper_model(Gpar,startsim, tprm, prm, ptqfile, utfile, r2fil,
        modstate, savestate,kal, spinup) # a single run 
    println(catchment)
    println("KGE=",round(KGE,digits=3))
    println("NSE=",round(NSE,digits=3))
    println("bias=",round(bias,digits=3))
end
if(kal == 1) # calibrate
    #                   u,        pro,         TX,        Pkorr,    skorr,          GscInt,         OVP     OVIP 
    param_range = [(1.0,3.0), (0.05,0.05), (-0.5, 0.5), (0.5, 2.0), (0.5,2.0), (0.065,0.075), (tprm[7],tprm[7]),
        (tprm[8],tprm[8]), (tprm[9],tprm[9]),(tprm[10],tprm[10])] # 
    println(param_range)
    calib_single_wsh_tmp(param) = calib_single_wsh(Gpar,startsim, param, prm, ptqfile, utfile, r2fil,
                                           modstate, savestate, kal, spinup)
    res = bboptimize(calib_single_wsh_tmp; SearchRange = param_range, MaxSteps = 1000, TraceMode = :verbose, NThreads=Threads.nthreads()-1)
    param_hydro = best_candidate(res)
    println(param_hydro)
end
println("Pkorr = ", round(tprm[4], digits=3))
println("Time elapsed = ", time() - t1, " seconds")
