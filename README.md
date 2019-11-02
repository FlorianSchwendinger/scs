# R interface to the [Splitting Conic Solver (SCS)](https://github.com/cvxgrp/scs)

## Control Arguments
| Parameter               | Type    |             Description                                      | Default Value |
|:------------------------|:--------|:-------------------------------------------------------------|:-------------:|
| `max_iters`             | integer | maximum number of iterations                                 |   `5000`      |
| `normalize`             | boolean | turn heuristic data rescaling on and off                     |   `TRUE`      |
| `verbose`               | boolean | write out progress                                           |   `FALSE`     |
| `cg_rate`               | numeric | for indirect, tolerance goes down like `(1 / iter)^cg_rate`  |      `2`      |
| `scale`                 | numeric | if normalized, rescales by this factor                       |      `1`      |
| `rho_x`                 | numeric | x equality constraint scaling                                |   `1e-3`      |
| `alpha`                 | numeric | relaxation parameter                                         |    `1.5`      |
| `eps`                   | numeric | convergence tolerance                                        |   `1e-5`      |
| `acceleration_lookback` | integer | number of iterations to look back for Anderson acceleration  |   `20`        |

## Reference
* O'Donoghue Brendan, Chu Eric, Parikh Neal, Boyd Stephen (2016).
  Operator Splitting for Conic Optimization via Homogeneous Self-Dual Embedding.
  Journal of Optimization Theory and Applications 169(3), 1042-1068.
  URL [`https://doi.org/10.1007/s10957-016-0892-3`](https://doi.org/10.1007/s10957-016-0892-3).
