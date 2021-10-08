## Summary of modifications

These modifications are against scs 3.0 (3aaa93c7aa04c7001df5e51b81f21b126dfa99b3).

1. In `rw.c`, all the `fread` calls were checked for returned `bytes`
   to avoid compilation warnings. Routines affected are
   `read_scs_cone`, `read_scs_stgs`, `read_amatrix`, `read_scs_data`,
   `SCS`. 

2. In `scs.mk` added a solaris check to ensure it compiles using
   Oracle developer studio or GNU tools. This is for CRAN.
   
   
