include("calibration.jl")
@assert length(ARGS) == 1
calibrateMultipleCatchments(ARGS[1])
