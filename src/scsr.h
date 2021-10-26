
#ifdef _OPENMP
    #include <omp.h>
#else
    #define omp_get_num_threads() -1
#endif
