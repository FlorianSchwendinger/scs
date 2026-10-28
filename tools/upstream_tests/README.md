# Regenerating the saved test problems

`inst/tinytest/problem_results/*.RDS` are recordings of the C test suite
that ships with the bundled SCS sources (`src/scs/test/problems/*.h`),
one file per `scs()` / `scs_solve()` call. `test_problems.R` replays them
through the package. This directory holds the tooling that produces them;
it is not part of the built package (`.Rbuildignore`).

    bash tools/upstream_tests/regenerate.sh

| file | role |
|---|---|
| `regenerate.sh` | the whole pipeline, run from the package root |
| `dump_tests.c` | includes upstream's test headers with `scs`, `scs_init`, `scs_solve` macro-redirected to wrappers that record data, cone, settings, warm start, and result, and `gen_random_prob_data` redirected to record the planted optimum (`c'x*`, `b'y*`); built with `R CMD SHLIB`, entry `dump_tests_run(outdir)` |
| `run_dump.R` | loads the shared object inside R and runs the tests from the source tree root (their data paths are relative) |
| `dumps_to_rds.R` | one record to one RDS in the layout `test_problems.R` reads |

Design points worth keeping:

- The library is built from a copy of `src/scs` with `src/scs_mods`
  applied and the flags of `src/Makevars`, linked against R's BLAS/LAPACK,
  and the tests run inside an R process, so the problem data are exactly
  what SCS saw upstream. The solution SCS returned is not stored: iterates
  differ across platforms, so `test_problems.R` checks the replayed
  solution with upstream's own assertions. When a file is added, its
  upstream assertions go into the `expected` table in `test_problems.R`;
  the test fails for a file with no entry.
- `rob_gauss_cov_est` call 3 only writes a data file and is not recorded.
- Doubles cross from C to R as hexadecimal floats (`%a`). R's parser does
  not round-trip 17-digit decimals exactly, and a one-ulp change in the
  data changes the iterates.
- Skipped as not expressible through the one-shot R API: re-solves on a
  reused workspace (SCS keeps the Anderson acceleration memory and the
  row scaling between `scs_solve()` calls), validation cases that corrupt
  the matrix dimensions or leave the settings unset, and the upstream
  tests for `scs_update`, file I/O, the spectral cones, and the internal
  helper functions. The list of tests run is the `RUN(...)` block at the
  end of `dump_tests.c`.
