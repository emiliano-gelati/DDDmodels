# Calibration tools for DDD


## 1. Requirements

Julia version and libraries are defined in `Project.toml` and `Manifest.toml`:
See [Pkg documentation](https://pkgdocs.julialang.org/v1/environments/#Using-someone-else's-project)
for how to recreate the Julia environment to calibrate DDD.


## 2. How to set upo and run a calibration

To calibrate DDD on a set of catchments using 10 parallel threads, run

`julia --project=. --threads 11 calibrate_DDD.jl <path/to/settings.toml>`

where the last argument is the path to a TOML file (for example `settings/calibration_251120.toml`)
containing the following entries:

- `root_output`: absolute path to the folder under which all calibration output will be stored.
- `path_catchments_list`: path to the text file with the IDs of the catchments to be calibrated,
  where lines starting with `#` are ignored (example: `settings/catchments/catchments_251120.txt`).
- `path_parameter_ranges`: path to the TOML file defining upper and lower bounds for the hydrologic
  parameters to be calibrated, where the `[default]` section is compulsory; and sections for
  individual catchments can be added by using catchment IDs for headers, where custom ranges for
  individual parameters can be set
  (example: `settings/parameter_ranges/parameter_ranges_251120.toml`).
- `template_path_ptq`: template path to the PTQ (precipitation, temperature and discharge) CSV file
  for each catchment, where `<CATCHMENT>` and `<PERIOD>` are replaced by the catchment ID and
  the abbreviation for each period (e.g. `kal` for calibration) in the function `pathsPTQ`.
- `template_path_inipar`:template path to the files with initial values for the full parameter set
  defining each catchment (example: files in `/hdata/hmdata01/DDD_calibration/BestPrimo2025/`).
- `spinup`: model spinup time (days).
- `steps_max`: maximum number of steps used by the optimisation algorithm.
- `periods`: section defining names (arbitrary) and abbreviations (must match the tokens used to
  replace `<PERIOD>` in `template_path_ptq`) for each period the model is run on.


## 3. Output

For each catchment, results are written in folder `root_output/<CATCHMENT>/`, where the `initial/`
and `calibrated` subfolders contains output from DDD runs with initial and calibrated parameter
values, respectively.

In each subfolder, the following files are written by
`../DDDFunctions/DDDAllTerrain22012024.jl`:

- `r2_<PERIOD>.csv`: row file with NSE, KGE, bias and hydrologic parameter values.
- `series_<PERIOD>.csv`: DDD output time series.

Additional files in the `calibrated/` subfolder:

- `parameters.csv`: new full parameter set defining the catchment, including the newly calibrated
  hydrologic parameters.
- `logged_solutions.csv`: KGE and hydrologic parameter values at the last iteration of the
  optimisation algorithm. This log is produced by
  [BlackBoxOptim.jl](https://docs.sciml.ai/Optimization/stable/optimization_packages/blackboxoptim/).


## 4. Possible improvements

1. Redefine DDDAll... as made up of 2 functions: f1 loading PTQ time series, and f2 doing the rest
   taking PTQ series among its input. To calibrate a catchment, f1 would be called only once and
   f2 would be iterated over by the search algorith.
2. Specify which parameters should be calibrated and which have a fixed value in settings,
   instead of using collapsed ranges, to reduce the formal number of parameters.
3. Check that functions are type-stable to avoid performance losses.
4. Profile DDD code to look for performance bottlenecks.
5. Test if using views instead of slices improves performance by reducing memory allocations.
6. DDDAllTerrain: select columns of ptqinn with syntax [:,"q"] so that the output is Vector?
7. Sensitivity analysis (e.g. [GlobalSensitivity.jl](https://docs.sciml.ai/GlobalSensitivity/stable/))
