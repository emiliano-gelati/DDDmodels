using Configurations
using Infiltrator
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
    path_catchments_list::String
    path_parameter_ranges::String
    template_path_ptq::String
    template_path_inipar::String
    spinup::Int
    steps_max::Int
    periods::Dict{String,String}
end

function pathsPTQ(id::String, settings::SettingsCalibration)
    Dict(k => replace(settings.template_path_ptq, "<CATCHMENT>" => id, "<PERIOD>" => v) for (k, v) in settings.periods)
end

function pathIniPar(id::String, settings::SettingsCalibration)
    replace(settings.template_path_inipar, "<CATCHMENT>" => id)
end

function checkInput(id::String, settings::SettingsCalibration)
    paths = pathsPTQ(id, settings)
    paths["inipar"] = pathIniPar(id, settings)
    Dict(vcat(["id" => id],[k => isfile(p) for (k, p) in paths]))
end

function runDDD(paths_ptq::Dict{String,String}, params_hyd::Vector{Float64}, params_all::DataFrame, spinup::Int, dir_out::String, txt::String)
    print("\tRuns with ", txt)
    kge_score = Dict(k => NaN for k in eachindex(paths_ptq))
    Threads.@threads for p in collect(eachindex(paths_ptq))
        path_out_series = joinpath(dir_out, "series_$(p).csv")
        path_out_r2 = joinpath(dir_out, "r2_$(p).csv")
        kge_score[p] = DDDAllTerrain(nothing, 1, params_hyd, params_all, paths_ptq[p], path_out_series, path_out_r2, 0, 0, 0, spinup, true)[3]
    end
    println(" -> KGE (period): ", join(["$(round(v, digits=4)) ($k)" for (k, v) in kge_score], ", "))
end

function makeEvaluator(parameters_all::DataFrame, path_ptq::String, spinup::Int, path_r2::String)
    function wrapper(hydpar::Vector{Float64})
        kge_score = DDDAllTerrain(nothing, 1, hydpar, parameters_all, path_ptq, "", path_r2, 0, 0, 1, spinup, true)[3]
        return 1. - kge_score
    end
end

function main(path_toml::String)
    # Names and positions (within the full set) of parameters to be calibrated 
    names_hydpar = ["u", "pro", "TX", "pkorr", "skorr", "GscInt", "OFVP", "OFVIP", "Lv", "Rv"]
    positions_hydpar = [20, 21, 22, 18, 19, 33, 34, 35, 36, 37]
    # Load settings from TOML file
    settings = from_toml(SettingsCalibration, path_toml)
    num_threads = max(Threads.nthreads() - 1, 1)
    # Load catchment list
    catchments = readlines(settings.path_catchments_list)
    filter!(line -> !isempty(strip(line)) && !startswith(strip(line), "#"), catchments)
    # Check that input files are valid
    ok_input = DataFrame(checkInput.(catchments, Ref(settings)))
    is_bad = (!all).(eachrow(ok_input[:,setdiff(names(ok_input), ["id"])]))
    if any(is_bad)
        error(println("Some input files do not exist or are not valid:\n", ok_input[is_bad,:]))
    end
    # Load parameter ranges
    raw = TOML.parsefile(settings.path_parameter_ranges)
    bounds_all = Dict(k => DataFrame(convert(Dict{String,Vector{Float64}}, d)) for (k, d) in raw)
    # Loop through catchments
    for (n, id) in enumerate(catchments)
        ## Root folder for catchment output
        dir_out = mkpath(joinpath(settings.root_output, id))
        println("\nCATCHMENT ", id, " (", n, " of ", length(catchments), "): output in ", dir_out)
        ## Paths to PTQ input for each period
        paths_ptq = pathsPTQ(id, settings)
        ## Load initial parameters and run DDD
        path_inipar = replace(settings.template_path_inipar, "<CATCHMENT>" => id)
        parameters_all = CSV.read(path_inipar, DataFrame, header=["Name", "val"], delim=';')
        parameters_hyd = [parameters_all.val[i] for i in positions_hydpar]
        dir_out_ini = mkpath(joinpath(dir_out, "initial"))
        runDDD(paths_ptq, parameters_hyd, parameters_all, settings.spinup, dir_out_ini, "initial parameters")
        ## Modify parameter bounds if specified for catchment (NOT IMPLEMENTED YET)
        if haskey(bounds_all, id)
            error("Not implemented yet: parameter ranges for individual catchments")
        else
            bounds_local = [Tuple(bounds_all["default"][:,k]) for k in names_hydpar]
        end
        ## Calibrate
        dir_out_cal = mkpath(joinpath(dir_out, "calibrated"))
        dir_log_cal = mkpath(joinpath(dir_out_cal, "log"))
        template_path_r2 = joinpath(dir_log_cal, "r2.csv")
        evaluator = makeEvaluator(parameters_all, paths_ptq["calibration"], settings.spinup, template_path_r2)
        print("\tCalibration started on ", Dates.format(now(), "yyyy-mm-dd HH:MM"))
        res = redirect_stdio(stdout=devnull, stderr=devnull) do
            bboptimize(evaluator; SearchRange=bounds_local, MaxSteps=settings.steps_max, TraceMode=:silent, SaveTrace=true, NThreads=num_threads)
        end
        print(" and ended after ", canonicalize(Second(Int(round(res.elapsed_time)))))
        println(" -> KGE: ", round(1 - best_fitness(res), digits=4))
        ## Write calibrated parameters to file
        parameters_hyd = best_candidate(res)
        parameters_all[positions_hydpar,"val"] .= parameters_hyd
        CSV.write(joinpath(dir_out_cal, "parameters.csv"), parameters_all, delim=';', writeheader=false)
        ## Merge calibration log files (1 r2fil per thread): TO DO (function in DDDAll... to rename r2fil and reuse it here)!
        ## Run DDD with calibrated parameters
        runDDD(paths_ptq, parameters_hyd, parameters_all, settings.spinup, dir_out_cal, "calibrated parameters")
    end
end

main(ARGS[1])
