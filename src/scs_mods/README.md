## Summary of modifications

These modifications are against [scs 3.3.1](https://github.com/cvxgrp/scs/tree/3.3.1)
(commit `cc244ebc`), the tag `src/scs` points to. `src/Makevars` and
`src/Makevars.win` copy the files in this directory over the upstream tree
before `libscsdir.a` is built; the submodule itself is never committed in
modified form. `tools/make_scs_mods.R` regenerates both files from the pristine
submodule and stops if upstream changed the text it edits, so run it after
every submodule bump.

1. `scs.mk`: the `OPT` default (`-O3 -march=native -fno-math-errno`) is
   emptied so that R's own `CFLAGS` govern (`-march=native` is not portable),
   and `-Wwrite-strings`, `-Wstrict-prototypes` and
   `-Werror=incompatible-pointer-types` are removed from the `CFLAGS` line
   (compiler specific; CRAN does not allow `-Werror`).

2. `cones.c`: the 12 `sprintf` calls in `format_cone_header()` become
   `snprintf` with the buffer size (`R CMD check` flags `sprintf` in
   compiled code).

3. `cones.c`: the complex PSD cone calls LAPACK `zheevr`, which R's own
   LAPACK does not provide (R ships `zheev` and `zheevd`), and "Writing R
   Extensions" allows a package to rely only on routines in R's LAPACK.
   The prototype, the workspace query and the projection call are switched
   to `zheev`; `zheev` returns the eigenvectors in place, so one `memcpy`
   moves them from `cXs` to `cZ` where the following code expects them.
   Without this the shared object fails to load on every R built with the
   internal LAPACK, for example CRAN's Windows machines. The real PSD cone
   keeps `dsyevr`, which R provides.
