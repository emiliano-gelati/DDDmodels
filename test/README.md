# Tests

## Basic test checking that model output is as expected

The test, which is run by executing `julia --project=.. runtests.jl`, checks whether,
given input parameters and PTQ time series, the model reproduces the expected output.
Data used in this test are in the `data/` subfolder:

- `parameters_12.70.csv`: model parameters
- `ptq_12.70_2009-2024.csv`: input PTQ time series
- `benchmark_series_12.70.csv`: expected output time series
- `output`: model output produced by the test (TO DO: automatic deletion upon successful test)

