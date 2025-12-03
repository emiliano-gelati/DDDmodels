# Calibration tools for DDD

## 1. Requirements

Julia version and libraries are defined in `Project.toml` and `Manifest.toml`:
See [Pkg documentation](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project)
for how to recreate the Julia environment to calibrate DDD.


## 2. 

To calibrate DDD on a catchment using 10 parallel threads, run in this folder

`julia --project=. --threads 11 calibrate_DDD.jl <path/to/settings.toml>`

where the last argument is the path to a TOML file (for example `settings/calibration_251120.toml`)
containing the following entries:

- `root_output`: absolute path to the folder under which all calibration output will be stored.
- `template_path_ptq`:
- `path_catchments_list`:
- `path_parameter_ranges`:
- `path_initial_parameters`:
- `spinup`: model spinup time (days).
- `steps_max`: maximum number of steps used by the optimisation algorithm.
- `periods`:

For each catchment, the following files are produced:

- ...


## Possible improvements

1. Redefine DDDAll... as made up of 2 functions: f1 loading PTQ time series, and f2 doing the rest
   taking PTQ series among its input. To calibrate a catchment, f1 would be called only once and
   f2 would be iterated over by the search algorith.
2. Specify which parameters should be calibrated and which have a fixed value in settings,
   instead of using collapsed ranges, to reduce the formal number of parameters?
3. Profile DDD code
4. Check that functions are type-stable
5. Test if using views instead of slices improves performance by reducing memory allocations
6. DDDAllTerrain: select columns of ptqinn with syntax [:,"q"] so that the output is Vector?
7. Sensitivity analysis (e.g. [GlobalSensitivity.jl](https://docs.sciml.ai/GlobalSensitivity/stable/))
