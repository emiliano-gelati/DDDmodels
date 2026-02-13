using Configurations
using Infiltrator
using TOML
using CSV
using Random
using Dates
using DataFrames
using BlackBoxOptim
include(joinpath(dirname(@__DIR__), "DDDFunctions", "DDDAllTerrain22012024.jl"))

Random.seed!(0)

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

mutable struct ParameterSet
    positions_hyd::Vector{UInt8}
    all::DataFrame
    hydrologic::Vector{Float64}
    function ParameterSet(path_in::String)
        positions_hyd::Vector{UInt8} = [20, 21, 22, 18, 19, 33, 34, 35, 36, 37]
        params_all = CSV.read(path_in, DataFrame, header=["Name", "val"], delim=';')
        hydrologic::Vector{Float64} = [params_all[i,"val"] for i in positions_hyd]
        new(positions_hyd, params_all, hydrologic)
    end
end

function setHydrologicParameters!(parameters::ParameterSet, hydrologic::Vector{Float64})
    parameters.all[parameters.positions_hyd,"val"] .= hydrologic
end

function pathsPTQ(id::String, settings::SettingsCalibration)
    Dict(k => replace(settings.template_path_ptq, "<CATCHMENT>" => id, "<PERIOD>" => v) for (k, v) in settings.periods)
end

function pathIniPar(id::String, settings::SettingsCalibration)
    replace(settings.template_path_inipar, "<CATCHMENT>" => id)
end

function dirCatchment(id::String, settings::SettingsCalibration)
    joinpath(settings.root_output, id)
end

function pathDone(id::String, settings::SettingsCalibration)
    joinpath(dirCatchment(id, settings), "done")
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
        kge_score[p] = DDDAllTerrain(fill(NaN, 2), 1, params_hyd, params_all, paths_ptq[p], path_out_series, path_out_r2, 0, 0, 0, spinup, true)[3]
    end
    println(" -> KGE (period): ", join(["$(round(v, digits=4)) ($k)" for (k, v) in kge_score], ", "))
end

function makeEvaluator(parameters_all::DataFrame, path_ptq::String, spinup::Int, path_r2::String)
    function wrapper(hydpar::Vector{Float64})
        kge_score = DDDAllTerrain(fill(NaN, 2), 1, hydpar, parameters_all, path_ptq, "", path_r2, 0, 0, 1, spinup, true)[3]
        return 1. - kge_score
    end
end

function calibrateMultipleCatchments(path_toml::String)
    # Load settings from TOML file
    settings = from_toml(SettingsCalibration, path_toml)
    num_threads = max(Threads.nthreads() - 1, 1)
    # Load catchment list and keep only those not done yet
    catchments = readlines(settings.path_catchments_list)
    filter!(line -> !isempty(strip(line)) && !startswith(strip(line), "#"), catchments)
    num_tot = length(catchments)
    filter!(id -> !isfile(pathDone(id, settings)), catchments)
    if num_tot > length(catchments)
        println("Of ", num_tot, " catchments, ", num_tot - length(catchments), " are already calibrated")
    end
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
        dir_out = mkpath(dirCatchment(id, settings))
        println("\nCATCHMENT ", id, " (", n, " of ", length(catchments), "): output in ", dir_out)
        ## Paths to PTQ input for each period
        paths_ptq = pathsPTQ(id, settings)
        ## Load initial parameters and run DDD
        path_inipar = replace(settings.template_path_inipar, "<CATCHMENT>" => id)
        parameters = ParameterSet(path_inipar)
        dir_out_ini = mkpath(joinpath(dir_out, "initial"))
        runDDD(paths_ptq, parameters.params_hyd, parameters.params_all, settings.spinup, dir_out_ini, "initial parameters")
        ## Modify parameter bounds if specified for catchment (NOT IMPLEMENTED YET)
        if haskey(bounds_all, id)
            error("Not implemented yet: parameter ranges for individual catchments")
        else
            names_hydpar = ["u", "pro", "TX", "pkorr", "skorr", "GscInt", "OFVP", "OFVIP", "Lv", "Rv"]
            bounds_local = [Tuple{Float64,Float64}(bounds_all["default"][:,k]) for k in names_hydpar]
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
        setHydrologicParameters(parameters, best_candidate(res))
        CSV.write(joinpath(dir_out_cal, "parameters.csv"), parameters.all, delim=';', writeheader=false)
        ## Merge calibration log files (1 r2fil per thread): TO DO (function in DDDAll... to rename r2fil and reuse it here)!
        ## Run DDD with calibrated parameters
        runDDD(paths_ptq, parameters.hydrologic, parameters.all, settings.spinup, dir_out_cal, "calibrated parameters")
        ## Create empty file ("done") to be used in case of restart to skip this catchment
        touch(pathDone(id, settings))
    end
end

function runSingleCatchment(path_toml::String, id::String, period::String)
    # Load settings from TOML file
    settings = from_toml(SettingsCalibration, path_toml)
    # Path to PTQ input
    path_ptq = pathsPTQ(id, settings)[period]
    # Load initial parameters
    path_inipar = replace(settings.template_path_inipar, "<CATCHMENT>" => id)
    parameters = ParameterSet(path_inipar)
    # Root folder for catchment output
    dir_out = mkpath(joinpath(settings.root_output, "single_runs", id))
    path_out_series = joinpath(dir_out, "series_$(id)_$(period).csv")
    path_out_r2 = joinpath(dir_out, "r2_$(id)_$(period).csv")
    println("Output in ", dir_out)
    # Run DDD
    DDDAllTerrain(fill(NaN, 2), 1, parameters.hydrologic, parameters.all, path_ptq, path_out_series, path_out_r2, 0, 0, 0, settings.spinup, true)
end
