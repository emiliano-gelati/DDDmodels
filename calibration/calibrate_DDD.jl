using Configurations
using TOML
using CSV
using DataFrames
using BlackBoxOptim
using Distributions
using LsqFit
using Statistics
using Dates
using JLD2
# Folder with DDD main program and functions
dir_DDD = joinpath(dirname(@__DIR__), "DDDFunctions")
# Preprocessing routines
include(joinpath(dir_DDD, "Big2SmallLambda.jl"))
include(joinpath(dir_DDD, "CeleritySubSurface.jl"))
include(joinpath(dir_DDD, "SingleUH.jl"))
include(joinpath(dir_DDD, "SingleNormalUH.jl"))
include(joinpath(dir_DDD, "LayerEstimation.jl"))
include(joinpath(dir_DDD, "PyrAreas.jl"))
include(joinpath(dir_DDD, "GrWPoint.jl"))
include(joinpath(dir_DDD, "RiverPoint.jl"))
include(joinpath(dir_DDD, "TemperatureVector.jl"))
# EB and Snow Routines
include(joinpath(dir_DDD, "NedbEBGlac_debug04072022.jl"))
include(joinpath(dir_DDD, "SnowpackTemp.jl"))
include(joinpath(dir_DDD, "TempstartUpdate.jl"))
include(joinpath(dir_DDD, "SmeltEBGlac_debug04072022.jl"))
include(joinpath(dir_DDD, "CloudCoverGlac_debug04072022.jl"))
include(joinpath(dir_DDD, "TssDewpoint.jl"))
include(joinpath(dir_DDD, "SolradTransAlbedoper_hrs_debug04072022.jl"))
include(joinpath(dir_DDD, "LongWaveRad_debug04072022.jl"))
include(joinpath(dir_DDD, "SensibleLatHeat_debug04072022.jl"))
include(joinpath(dir_DDD, "AlbedoUEB_debug04072022.jl"))
include(joinpath(dir_DDD, "GroundPrecCC.jl"))
include(joinpath(dir_DDD, "SnowGamma.jl"))
include(joinpath(dir_DDD, "Varc.jl"))
include(joinpath(dir_DDD, "NewSnowDensityEB.jl"))
include(joinpath(dir_DDD, "NewSnowSDEB.jl"))
include(joinpath(dir_DDD, "DensityAge.jl"))
# Subsurface and Evaporation routines
include(joinpath(dir_DDD, "LayerCapacityUpdate.jl"))
include(joinpath(dir_DDD, "PotentialEvapPT.jl"))
include(joinpath(dir_DDD, "UnsaturatedEvapEB.jl"))
include(joinpath(dir_DDD, "LayerEvap.jl"))
include(joinpath(dir_DDD, "UnsaturatedExEvap.jl"))
include(joinpath(dir_DDD, "WetlandsEB.jl"))
include(joinpath(dir_DDD, "GrvInputDistributionICap2022.jl"))
include(joinpath(dir_DDD, "OFICap.jl"))
include(joinpath(dir_DDD, "LayerUpdate.jl"))
include(joinpath(dir_DDD, "BogLayerUpdate.jl"))
include(joinpath(dir_DDD, "RiverUpdate.jl"))
# Overland Flow routine
include(joinpath(dir_DDD, "OverlandFlowDynamicDD.jl"))
# Efficiency criteria
include(joinpath(dir_DDD, "NSE_ths.jl"))
include(joinpath(dir_DDD, "KGE_ths.jl"))
# DDD main program
include(joinpath(dir_DDD, "DDDAllTerrain22012024.jl"))

@option struct SettingsCalibration
    root_output::String
    template_path_ptq::String
    path_catchments_list::String
    path_parameter_ranges::String
    path_initial_parameters::String
    spinup::Int
    steps_max::Int
    periods::Dict{String,String}
end

function runDDD(paths_ptq::Dict{String,String}, params_hyd::Vector{Float64}, params_all::DataFrame, spinup::Int, dir_out::String)
    kge_score = Dict(k => NaN for k in eachindex(paths_ptq))
    runtime = Dict(k => NaN for k in eachindex(paths_ptq))
    Threads.@threads for p in collect(eachindex(paths_ptq))
        path_out_series = joinpath(dir_out, "series_$(p).csv")
        path_out_r2 = joinpath(dir_out, "r2_$(p).csv")
        t0 = time()
        kge_score[p] = DDDAllTerrain(nothing, 1, params_hyd, params_all, paths_ptq[p], path_out_series, path_out_r2, 0, 0, 0, spinup, true)[3]
        runtime[p] = time() - t0
    end
    println("    - KGE: ", join(["$k ($(round(v, digits=4)))" for (k, v) in kge_score], ", "))
    println("    - Run times (s): ", join(["$k ($(round(Int, v)))" for (k, v) in runtime], ", "))
end

function makeEvaluator(parameters_all::DataFrame, path_ptq::String, spinup::Int)
    function wrapper(hydpar::Vector{Float64})
        kge_score = DDDAllTerrain(nothing, 1, hydpar, parameters_all, path_ptq, "", "", 0, 0, 1, spinup, true)[3]
        return 1 - kge_score
    end
    return wrapper
end

function main(path_toml::String)
    # Names and positions of (hydrologic) parameters to be calibrated in the full parameter set
    names_hydpar = ["u", "pro", "TX", "pkorr", "skorr", "GscInt", "OFVP", "OFVIP", "Lv", "Rv"]
    positions_hydpar = [20, 21, 22, 18, 19, 33, 34, 35, 36, 37]
    # Load settings from TOML file
    settings = from_toml(SettingsCalibration, path_toml)
    # Load catchment list
    catchments = readlines(settings.path_catchments_list)
    filter!(line -> !isempty(strip(line)) && !startswith(strip(line), "#"), catchments)
    # Load parameter ranges
    raw = TOML.parsefile(settings.path_parameter_ranges)
    bounds_all = Dict(k => DataFrame(convert(Dict{String,Vector{Float64}}, d)) for (k, d) in raw)
    # Loop through catchments
    for (n, id) in enumerate(catchments)
        println("\nCATCHMENT ", id, "(", n, "/", length(catchments), ")")
        ## Paths to PTQ input for each period
        paths_ptq = Dict(k => replace(settings.template_path_ptq, "<CATCHMENT>" => id, "<PERIOD>" => v) for (k, v) in settings.periods)
        ## Root folder for catchment output
        dir_out = mkpath(joinpath(settings.root_output, id))
        ## Load initial parameters and run DDD
        path_inipar = replace(settings.path_initial_parameters, "<CATCHMENT>" => id)
        parameters_all = CSV.read(path_inipar, DataFrame, header=["Name", "val"], delim=';')
        parameters_hyd = [parameters_all.val[i] for i in positions_hydpar]
        dir_out_ini = mkpath(joinpath(dir_out, "initial"))
        println("  * DDD runs using initial parameters: output in ", dir_out_ini)
        runDDD(paths_ptq, parameters_hyd, parameters_all, settings.spinup, dir_out_ini)
        ## Calibrate
        num_threads = max(Threads.nthreads() - 1, 1)
        if haskey(bounds_all, id)
            error("Not implemented yet: parameter ranges for individual catchments")
        else
            bounds_local = [Tuple(bounds_all["default"][:,k]) for k in names_hydpar]
        end
        evaluator = makeEvaluator(parameters_all, paths_ptq["calibration"], settings.spinup)
        println("  * Calibration started on ", now())
        res = bboptimize(evaluator; SearchRange=bounds_local, MaxSteps=settings.steps_max, TraceMode=:silent, NThreads=num_threads)
        println("  * Calibration ended after ", Dates.canonicalize(Dates.CompoundPeriod(Second(Int(round(res.elapsed_time))))))
        println("  * Best KGE: ", round(1 - best_fitness(res), digits=4))
        ## Write calibrated parameters to file
        parameters_hyd = best_candidate(res)
        parameters_all[positions_hydpar,"val"] .= parameters_hyd
        dir_out_cal = mkpath(joinpath(dir_out, "calibrated"))
        path_out_calpar = joinpath(dir_out_cal, "parameters_$(id).csv")
        CSV.write(path_out_calpar, parameters_all, delim=';', writeheader=false)
        println("  * Calibrated parameters saved to ", path_out_calpar)
        ## TO DO: Save whole calibration history (parameter values and KGE at each iteration) to file
        ## Run DDD with calibrated parameters
        println("  * DDD runs using calibrated parameters: output in ", dir_out_cal)
        runDDD(paths_ptq, parameters_hyd, parameters_all, settings.spinup, dir_out_cal)
    end
end

main(ARGS[1])

# TO CHECK:
# - DDDAllTerrain: select columns of ptqinn with syntax [:,"q"] so that the output is Vector?
