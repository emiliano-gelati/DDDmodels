# Calibration tools for DDD

To calibrate DDD on a catchment using 10 parallel threads, run in this folder

`julia --project=. --threads 11 calibrate_DDD.jl settings/calibration_251120.toml`

where the last argument is the path to a TOML file ...

Julia version and dependencies are defined in `Project.toml` and `Manifest.toml`:
See [Pkg documentation](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project)
for how to recreate the Julia environment to calibrate DDD.
