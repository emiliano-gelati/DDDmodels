# To run this script:
# julia --project=. run_single_catchment.jl settings/calibration_251120.toml 2.11 calibration
include("calibration.jl")

@assert length(ARGS) == 3
path_settings, id, period = ARGS

# Simple run
#runSingleCatchment(path_settings, id, period)

#using Profile
#Profile.clear()
#@profile runSingleCatchment(path_settings, id, period)
#Profile.print(format=:tree, sortedby=:count, maxdepth=10, mincount=5)

#using InteractiveUtils
#@code_warntype runSingleCatchment(path_settings, id, period)

using BenchmarkTools
@benchmark runSingleCatchment(path_settings, id, period)

# ProfileView displays works in REPL and not script
#include("calibration.jl")
#using ProfileView
#using Cthulhu
#@profview runSingleCatchment("settings/calibration_251120.toml", "2.11", "calibration")
