## Summary of modifications

These modifications are against [scs@fa32d87](https://github.com/cvxgrp/scs). 

1. In `rw.c`, all the `fread` calls were checked for returned `bytes`
   to avoid compilation warnings. Routines affected are
   `read_scs_cone`, `read_scs_stgs`, `read_amatrix`, `read_scs_data`,
   `SCS`. This was done by adding two macros `FREAD` and `FWRITE` and
   globally replacing `fread` with `FREAD` and `fwrite` with
   `FWRITE`, except for line 144.

2. Remove all SOLARIS patches as the platform is no longer supported.    
