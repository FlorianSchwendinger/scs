# R interface to the [Splitting Conic Solver (SCS)](https://github.com/cvxgrp/scs)

## Control Arguments


|Parameter             |Default Value |Type    |Description                                                                    |
|:---------------------|:-------------|:-------|:------------------------------------------------------------------------------|
|max_iters             |100000        |integer |maximum number of iterations                                                   |
|eps_rel               |1e-04         |numeric |relative feasibility tolerance                                                 |
|eps_abs               |1e-04         |numeric |absolute feasibility tolerance                                                 |
|eps_infeas            |1e-07         |numeric |Infeasibility tolerance (primal and dual)                                      |
|alpha                 |1.5           |numeric |Douglas-Rachford relaxation parameter                                          |
|rho_x                 |1e-06         |numeric |Primal scale factor                                                            |
|scale                 |0.1           |numeric |Initial dual scale factor, updated if `adaptive_scale` is True                 |
|verbose               |FALSE         |boolean |Write out progress                                                             |
|normalize             |TRUE          |boolean |Whether to perform heuristic data rescaling                                    |
|warm_start            |FALSE         |boolean |Whether to warm start                                                |
|acceleration_lookback |0             |integer |How much memory to use for Anderson acceleration                               |
|acceleration_interval |1             |integer |Run Anderson acceleration every this number of iterations                      |
|adaptive_scale        |TRUE          |boolean |Whether to heuristically adapt dual `scale` through the solve                  |
|log_csv_filename      |              |string  |If non `NULL` SCS will write logs of quantities through the solver (ignored)   |
|time_limit_secs       |0             |numeric |Time limit for solve run in seconds (can be fractional, 0 implies no limit     |
|write_data_filename   |              |string  |If non `NULL`set the problem data is dumped to this filename (ignored)         |


## Reference
* O'Donoghue Brendan, Chu Eric, Parikh Neal, Boyd Stephen (2016).
  Operator Splitting for Conic Optimization via Homogeneous Self-Dual Embedding.
  Journal of Optimization Theory and Applications 169(3), 1042-1068.
  URL [`https://doi.org/10.1007/s10957-016-0892-3`](https://doi.org/10.1007/s10957-016-0892-3).
