#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>
#include <Rinternals.h>

extern SEXP scsr(SEXP data, SEXP cone, SEXP params);
extern SEXP openmp_get_number_of_threads();

static const R_CallMethodDef CallEntries[] = {
					      {"scsr", (DL_FUNC) &scsr, 3},
					      {"openmp_get_number_of_threads", (DL_FUNC) &openmp_get_number_of_threads, 0},
					      {NULL, NULL, 0}
};

void R_init_scs(DllInfo *dll) {
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
