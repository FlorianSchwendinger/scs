# scs 3.3.1

- Synced up to version 3.3.1 of the SCS C library. The package version
  now matches the library version.
- New cones: complex positive semidefinite (`cs`). Box cone `bsize`
  defaults to `length(bl) + 1`.
- New `info` fields from 3.3.1, including `aa_stats`; new `scs_version()`.
- Anderson acceleration is now on by default, as in SCS:
  `acceleration_lookback = 10L`, `acceleration_interval = 5L` (were `0L`
  and `1L`). Pass the old values to restore the previous behavior.
- Unknown `control` names now give a warning (they were silently
  ignored); unknown cone fields are an error.
- A negative `acceleration_lookback` is now an error; use
  `acceleration_type_1 = FALSE` for type-II Anderson acceleration.
- `A` and `P` accept base matrices, `dgCMatrix` (and `dsCMatrix` for `P`)
  and `simple_triplet_matrix`; other `Matrix` classes give an error.
- Test suite regenerated from the upstream C tests (83 reference problems, checked with upstream's assertions).
- Upstream modifications are documented in `src/scs_mods/README.md`.

# scs 3.0-0

- Set `acceleration_lookback` parameter to (current default) `10L`
  from previous value of `20L`.
- Fixed up flags for compiling with Oracle Developer Studio
  compilers.
- Edited and cleaned up `test_socp` which had copy of `test_psd`
  repeated.

# scs 1.3-1

- Synced up to version scs-2.1.1 C library.
