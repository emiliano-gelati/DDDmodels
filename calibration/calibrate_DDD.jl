using Infiltrator
using TOML
using CSV
using Distributions
using LsqFit
using Statistics
using Dates
using DataFrames
using Plots
using BlackBoxOptim
using JLD2

path_dddfun = joinpath(dirname(@__DIR__), "DDDFunctions")
# Preprocessing routines
include(joinpath(path_dddfun, "Big2SmallLambda.jl"))
include(joinpath(path_dddfun, "CeleritySubSurface.jl"))
include(joinpath(path_dddfun, "SingleUH.jl"))
include(joinpath(path_dddfun, "SingleNormalUH.jl"))
include(joinpath(path_dddfun, "LayerEstimation.jl"))
include(joinpath(path_dddfun, "PyrAreas.jl"))
include(joinpath(path_dddfun, "GrWPoint.jl"))
include(joinpath(path_dddfun, "RiverPoint.jl"))
include(joinpath(path_dddfun, "TemperatureVector.jl"))
# EB and Snow Routines
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
# Subsurface and Evaporation routines
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
# Overland Flow routine
include(joinpath(path_dddfun, "OverlandFlowDynamicDD.jl"))
# Efficiency criteria
include(joinpath(path_dddfun, "NSE_ths.jl"))
include(joinpath(path_dddfun, "KGE_ths.jl"))
# Model Module
include(joinpath(path_dddfun, "DDDAllTerrain22012024.jl"))

#function makeEvaluator(parameters_initial::Vector{Float64}, path_ptq::String, spinup::Int)
#    return let
#        parameters_initial = parameters_initial,
#        path_ptq = path_ptq,
#        spinup = spinup
#        function wrapper(x::Vector{Float64})
#            return DDDAllTerrain(1, x, parameters_initial, path_ptq, "", "", 0, 0, 1, spinup)[3]
#        end
#    end
#end

function makeEvaluator(parameters_initial::DataFrame, path_ptq::String, spinup::Int)
    function wrapper(x::Vector{Float64})
        score_kge = DDDAllTerrain(1, x, parameters_initial, path_ptq, "", "", 0, 0, 1, spinup)[3]
        return 1 - score_kge
    end
    return wrapper
end


num_threads = max(Threads.nthreads() - 1, 1)
order_parameters = ["u", "pro", "TX", "Pkorr", "skorr", "GscInt", "OVP", "OVIP", "Lv", "rv"]
# Load settings from file
settings = TOML.parsefile(ARGS[1])
settings["periods"] = convert(Dict{String,String}, settings["periods"])
settings["path"] = convert(Dict{String,String}, settings["path"])
# Load catchment list
catchments = readlines(settings["path"]["catchments"])
filter!(line -> !isempty(strip(line)) && !startswith(strip(line), "#"), catchments)
# Load parameter ranges
raw = TOML.parsefile(settings["path"]["parameter_ranges"])
parameter_ranges = Dict(k => DataFrame(convert(Dict{String,Vector{Float64}}, d)) for (k, d) in raw)
# Loop through catchments
for id in catchments
    if haskey(parameter_ranges, id)
        error("Not implemented yet: parameter ranges for individual catchments")
    else
        bounds = [Tuple(parameter_ranges["default"][:,k]) for k in order_parameters]
    end
    path_ptq = replace(settings["path"]["ptq"], "{CATCHMENT}" => id, "{PERIOD}" => settings["periods"]["calibration"])
    parameters_initial = CSV.read(replace(settings["path"]["parameters_initial"], "{CATCHMENT}" => id),
                                  DataFrame, header=["Name", "val"], delim=';')
    evaluator = makeEvaluator(parameters_initial, path_ptq, settings["spinup"])
    res = bboptimize(evaluator; SearchRange=bounds, MaxSteps=1000, TraceMode=:verbose, NThreads=num_threads)
    param_hydro = best_candidate(res)
    dir_out = mkpath(joinpath(settings["path"]["output"], id))
    exit()
end
