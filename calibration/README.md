# Calibration tools for DDD

To calibrate DDD on a catchment using 10 parallel threads, run in this folder

`julia --project=. --threads 11 calibrate_DDD.jl`

Catchment selection and other settings are hardcoded in `calibrate_DDD.jl`,
which a modified version of `../DDDEcco/RunDDDv2.ipynb`, but these will be
set in a configuration file in the future.

Julia version and dependencies are defined in `Project.toml` and `Manifest.toml`:
See [Pks's documentation](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project)
for how to recreate the Julia environment to calibrate DDD.
