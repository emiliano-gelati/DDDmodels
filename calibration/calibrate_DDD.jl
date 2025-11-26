using Infiltrator
using TOML
using CSV
using DataFrames
using BlackBoxOptim
include(joinpath(dirname(@__DIR__), "DDDFunctions", "DDDAllTerrain22012024.jl"))

function makeEvaluator(parameters_initial::Vector{Float64}, path_ptq::String, spinup::Int)
    function wrapper(x::Vector{Float64})
        kge = DDDAllTerrain(1, x, parameters_initial, path_ptq, "", "", 0, 0, 1, spinup)[3]
        return 1 - kge
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
    path_initial_param = replace(settings["path"]["parameters_initial"], "{CATCHMENT}" => id)
    parameters_initial = CSV.read(path_initial_param, DataFrame, header=["Name", "val"], delim=';')
    # TESTS - START
    x = [parameters_initial.val[i] for i in [20, 21, 22, 18, 19, 33, 34, 35, 36, 37]]
    for i in range(1, 2)
        @time aux = DDDAllTerrain(1, x, parameters_initial.val, path_ptq, "", "", 0, 0, 1, settings["spinup"]);
    end
    @infiltrate
    # TESTS - END
    evaluator = makeEvaluator(parameters_initial.val, path_ptq, settings["spinup"])
    res = bboptimize(evaluator; SearchRange=bounds, MaxSteps=1000, TraceMode=:compact, NThreads=num_threads)
    param_hydro = best_candidate(res)
    dir_out = mkpath(joinpath(settings["path"]["output"], id))
    exit()
end
