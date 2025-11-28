using Infiltrator
using Configurations
using TOML
using CSV
using DataFrames
using BlackBoxOptim
include(joinpath(dirname(@__DIR__), "DDDFunctions", "DDDAllTerrain22012024.jl"))

@option struct SettingsCalibration
    name::String
    spinup::Int
    periods::Dict{String,String}
    path::Dict{String,String}
end

function runDDD(paths_ptq::Dict{String,String}, params_hyd::Vector{Float64}, params_all::Vector{Float64}, spinup::Int, dir_out::String)
    kge = Dict(k => NaN for k in eachindex(paths_ptq))
    runtimes = Dict(k => NaN for k in eachindex(paths_ptq))
    Threads.@threads for p in collect(eachindex(paths_ptq))
        path_out_series = joinpath(dir_out, "series_$(p).csv")
        path_out_r2 = joinpath(dir_out, "r2_$(p).csv")
        t0 = time()
        kge[p] = DDDAllTerrain(1, params_hyd, params_all, paths_ptq[p], path_out_series, path_out_r2, 0, 0, 0, spinup, true)[3]
        runtimes[p] = time() - t0
    end
    println("    - KGE: " * join(["$k ($(round(v, digits=4)))" for (k, v) in kge], ", "))
    println("    - Run times (s): " * join(["$k ($(round(Int, v)))" for (k, v) in runtimes], ", "))
end

function makeEvaluator(parameters_all::Vector{Float64}, path_ptq::String, spinup::Int)
    function wrapper(hydpar::Vector{Float64})
        kge = DDDAllTerrain(1, hydpar, parameters_all, path_ptq, "", "", 0, 0, 1, spinup, true)[3]
        return 1 - kge
    end
    return wrapper
end

function main(path_toml::String)
    # Names and positions of (hydrologic) parameters to be calibrated in the full parameter set
    names_hydpar = ["u", "pro", "TX", "Pkorr", "skorr", "GscInt", "OVP", "OVIP", "Lv", "rv"]
    positions_hydpar = [20, 21, 22, 18, 19, 33, 34, 35, 36, 37]
    # Load settings from TOML file
    settings = from_toml(SettingsCalibration, path_toml)
    # Load catchment list
    catchments = readlines(settings.path["catchments"])
    filter!(line -> !isempty(strip(line)) && !startswith(strip(line), "#"), catchments)
    # Load parameter ranges
    raw = TOML.parsefile(settings.path["parameter_ranges"])
    bounds_all = Dict(k => DataFrame(convert(Dict{String,Vector{Float64}}, d)) for (k, d) in raw)
    # Loop through catchments
    for (n, id) in enumerate(catchments)
        println("\nCATCHMENT $(id) ($(n)/$(length(catchments)))")
        ## Paths to PTQ input for each period
        paths_ptq = Dict(k => replace(settings.path["ptq"], "{CATCHMENT}" => id, "{PERIOD}" => v) for (k, v) in settings.periods)
        ## Root folder for catchment output
        root_out = mkpath(joinpath(settings.path["output"], id))
        ## Load initial parameters and run DDD
        path_inipar = replace(settings.path["parameters_initial"], "{CATCHMENT}" => id)
        parameters_all = CSV.read(path_inipar, DataFrame, header=["Name", "val"], delim=';')
        parameters_hyd = [parameters_all.val[i] for i in positions_hydpar]
        dir_out_ini = mkpath(joinpath(root_out, "initial"))
        println("  DDD runs using initial parameters (output in $(dir_out_ini)):")
        runDDD(paths_ptq, parameters_hyd, parameters_all.val, settings.spinup, dir_out_ini)
        ## Calibrate
        num_threads = max(Threads.nthreads() - 1, 1)
        if haskey(bounds_all, id)
            error("Not implemented yet: parameter ranges for individual catchments")
        else
            bounds_local = [Tuple(bounds_all["default"][:,k]) for k in names_hydpar]
        end
        evaluator = makeEvaluator(parameters_all.val, paths_ptq["calibration"], settings.spinup)
        res = bboptimize(evaluator; SearchRange=bounds_local, MaxSteps=1, TraceMode=:compact, NThreads=num_threads)
        ## Write calibrated parameters to file
        dir_out_cal = mkpath(joinpath(root_out, "calibrated"))
        parameters_hyd = best_candidate(res)
        findall(in(names_hydpar), parameters_all[:,"Name"]) # FIX: name discrepancies between all parameters and parameter ranges!
        ## Run DDD with calibrated parameters
        println("  DDD runs using calibrated parameters (output in $(dir_out_ini)):")
        runDDD(paths_ptq, parameters_hyd, parameters_all.val, settings.spinup, dir_out_ini)
        @infiltrate
        res.elapsed_time
    end
end

main(ARGS[1])

# TO DO:
# 1. select columns of ptqinn with syntax [:,"q"] so that the output is Vector
