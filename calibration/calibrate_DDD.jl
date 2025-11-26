using Infiltrator
using TOML
using CSV
using DataFrames
using BlackBoxOptim
include(joinpath(dirname(@__DIR__), "DDDFunctions", "DDDAllTerrain22012024.jl"))

function makeEvaluator(parameters_initial::Vector{Float64}, path_ptq::String, spinup::Int)
    function wrapper(x::Vector{Float64})
        kge = DDDAllTerrain(1, x, parameters_initial, path_ptq, "", "", 0, 0, 1, spinup, true)[3]
        return 1 - kge
    end
    return wrapper
end

function main(path_toml::String)
    # Load settings from TOML file
    settings = TOML.parsefile(path_toml)
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
        println("\nCATCHMENT $(id)")
        ## Output root folder for catchment
        dir_out = mkpath(joinpath(settings["path"]["output"], id))
        ## Load initial parameters
        path_initial_param = replace(settings["path"]["parameters_initial"], "{CATCHMENT}" => id)
        parameters_initial = CSV.read(path_initial_param, DataFrame, header=["Name", "val"], delim=';')
        x = [parameters_initial.val[i] for i in [20, 21, 22, 18, 19, 33, 34, 35, 36, 37]]
        ## Run DDD for all periods using initial parameters, and measure execution time
        dir_out_initial = mkpath(joinpath(dir_out, "initial"))
        println("Simulation results using initial parameters in $(dir_out_initial)")
        kge = Dict(k => NaN for k in eachindex(settings["periods"]))
        runtimes = Dict(k => NaN for k in eachindex(settings["periods"]))
        Threads.@threads for p in collect(eachindex(settings["periods"]))
            path_ptq = replace(settings["path"]["ptq"], "{CATCHMENT}" => id, "{PERIOD}" => settings["periods"][p])
            path_series = joinpath(dir_out_initial, "series_$(p)_$(id).csv")
            path_r2 = joinpath(dir_out_initial, "r2_$(p)_$(id).csv")
            t0 = time()
            kge[p] = DDDAllTerrain(1, x, parameters_initial.val, path_ptq, path_series, path_r2, 0, 0, 0, settings["spinup"], true)[3]
            runtimes[p] = time() - t0
        end
        println("KGE using initial parameters: " * join(["$k ($(round(v, digits=4)))" for (k, v) in kge], ", "))
        println("Run times (s): " * join(["$k ($(round(Int, v)))" for (k, v) in runtimes], ", "))
        ## Calibrate
        num_threads = max(Threads.nthreads() - 1, 1)
        order_parameters = ["u", "pro", "TX", "Pkorr", "skorr", "GscInt", "OVP", "OVIP", "Lv", "rv"]
        if haskey(parameter_ranges, id)
            error("Not implemented yet: parameter ranges for individual catchments")
        else
            bounds = [Tuple(parameter_ranges["default"][:,k]) for k in order_parameters]
        end
        path_ptq = replace(settings["path"]["ptq"], "{CATCHMENT}" => id, "{PERIOD}" => settings["periods"]["calibration"])
        evaluator = makeEvaluator(parameters_initial.val, path_ptq, settings["spinup"])
        res = bboptimize(evaluator; SearchRange=bounds, MaxSteps=1000, TraceMode=:compact, NThreads=num_threads)
        param_hydro = best_candidate(res)
        exit()
    end
end

main(ARGS[1])

# TO DO:
# 1. select columns of ptqinn with syntax [:,"q"] so that the output is Vector
