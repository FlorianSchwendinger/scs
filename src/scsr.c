#include <R.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include "scs/include/glbopts.h"
#include "scs/include/scs.h"
#include "scs/include/util.h"
#include "scs/linsys/amatrix.h"



SEXP c_to_r_integer(int c_int) {
    SEXP r_val;
    PROTECT(r_val = allocVector(INTSXP, 1));
    INTEGER(r_val)[0] = c_int;
    UNPROTECT(1);
    return r_val;
}

SEXP c_to_r_double(double c_double) {
    SEXP r_val;
    PROTECT(r_val = allocVector(REALSXP, 1));
    REAL(r_val)[0] = c_double;
    UNPROTECT(1);
    return r_val;
}

SEXP c_to_r_string(const char *c_string){
    SEXP r_val;
    PROTECT(r_val = allocVector(STRSXP, 1));
    SET_STRING_ELT(r_val, 0, mkChar(c_string));
    UNPROTECT(1);
    return r_val;
}

SEXP c_to_r_double_vec(double *c_double_vec, long len) {
    SEXP r_val;
    PROTECT(r_val = allocVector(REALSXP, len));
    for (long i = 0; i < len; i++) {
        REAL(r_val)[i] = c_double_vec[i];
    }
    UNPROTECT(1);
    return r_val;
}



SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++) {
    if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }
  return elmt;
}

scs_float getFloatFromListWithDefault(SEXP list, const char *str, scs_float def) {
    SEXP val = getListElement(list, str);
    if (val == R_NilValue) {
        return def;
    }
    val = coerceVector(val, REALSXP);
    return REAL(val)[0];
}

scs_int getIntFromListWithDefault(SEXP list, const char *str, scs_int def) {
    SEXP val = getListElement(list, str);
    if (val == R_NilValue) {
        return def;
    }
    val = coerceVector(val, INTSXP);
    return INTEGER(val)[0];
}

scs_float * getFloatVectorFromList(SEXP list, const char *str, scs_int * len) {
    SEXP vec = getListElement(list, str);
    *len = length(vec);
    vec = coerceVector(vec, REALSXP);
    return REAL(vec);
}

scs_int * getIntVectorFromList(SEXP list, const char *str, scs_int * len) {
    SEXP vec = getListElement(list, str);
    *len = length(vec);
    vec = coerceVector(vec, INTSXP);
    return INTEGER(vec);
}


SEXP info_to_r_list(ScsInfo * info) {
    SEXP r_info, r_info_names;

    PROTECT(r_info = NEW_LIST(12));
    PROTECT(r_info_names = NEW_CHARACTER(12));

    SET_STRING_ELT(r_info_names, 0, mkChar("iter"));
    SET_VECTOR_ELT(r_info, 0, c_to_r_integer(info->iter));

    SET_STRING_ELT(r_info_names, 1, mkChar("status"));
    SET_VECTOR_ELT(r_info, 1, c_to_r_string(info->status));

    SET_STRING_ELT(r_info_names, 2, mkChar("statusVal"));
    SET_VECTOR_ELT(r_info, 2, c_to_r_integer(info->status_val));

    SET_STRING_ELT(r_info_names, 3, mkChar("pobj"));
    SET_VECTOR_ELT(r_info, 3, c_to_r_double(info->pobj));

    SET_STRING_ELT(r_info_names, 4, mkChar("dobj"));
    SET_VECTOR_ELT(r_info, 4, c_to_r_double(info->dobj));

    SET_STRING_ELT(r_info_names, 5, mkChar("resPri"));
    SET_VECTOR_ELT(r_info, 5, c_to_r_double(info->res_pri));

    SET_STRING_ELT(r_info_names, 6, mkChar("resDual"));
    SET_VECTOR_ELT(r_info, 6, c_to_r_double(info->res_dual));

    SET_STRING_ELT(r_info_names, 7, mkChar("resInfeas"));
    SET_VECTOR_ELT(r_info, 7, c_to_r_double(info->res_infeas));

    SET_STRING_ELT(r_info_names, 8, mkChar("resUnbdd"));
    SET_VECTOR_ELT(r_info, 8, c_to_r_double(info->res_unbdd));

    SET_STRING_ELT(r_info_names, 9, mkChar("relGap"));
    SET_VECTOR_ELT(r_info, 9, c_to_r_double(info->rel_gap));

    SET_STRING_ELT(r_info_names, 10, mkChar("setupTime"));
    SET_VECTOR_ELT(r_info, 10, c_to_r_double(info->setup_time));

    SET_STRING_ELT(r_info_names, 11, mkChar("solveTime"));
    SET_VECTOR_ELT(r_info, 11, c_to_r_double(info->solve_time));

    SET_NAMES(r_info, r_info_names);
    UNPROTECT(2);

    return r_info;
}



SEXP scsr(SEXP n_variables, SEXP n_constraints, SEXP objective, 
          SEXP Ai, SEXP Ap, SEXP Ax, SEXP rhs, SEXP cone, SEXP params) {
    SEXP ret, retnames;

    /* allocate memory */
    ScsData * d = scs_malloc(sizeof(ScsData));
    ScsCone * k = scs_malloc(sizeof(ScsCone));
    ScsSettings * stgs = scs_malloc(sizeof(ScsSettings));
    ScsMatrix * A = scs_malloc(sizeof(ScsMatrix));
    ScsInfo * info = scs_calloc(1, sizeof(ScsInfo));
    ScsSolution * sol = scs_calloc(1, sizeof(ScsSolution));

    d->b = REAL(rhs);
    d->c = REAL(objective);
    d->n = INTEGER(n_variables)[0];
    d->m = INTEGER(n_constraints)[0];

    A->m = d->m;
    A->n = d->n;
    A->x = REAL(Ax);
    A->i = INTEGER(Ai);
    A->p = INTEGER(Ap);
    d->A = A;
  
    stgs->max_iters  = getIntFromListWithDefault(params, "max_iters", MAX_ITERS);
    stgs->normalize  = getIntFromListWithDefault(params, "normalize", NORMALIZE);
    stgs->verbose    = getIntFromListWithDefault(params, "verbose", VERBOSE);
    stgs->cg_rate    = getFloatFromListWithDefault(params, "cg_rate", CG_RATE);
    stgs->scale      = getFloatFromListWithDefault(params, "scale", SCALE);
    stgs->rho_x      = getFloatFromListWithDefault(params, "rho_x", RHO_X);
    stgs->alpha      = getFloatFromListWithDefault(params, "alpha", ALPHA);
    stgs->eps        = getFloatFromListWithDefault(params, "eps", EPS);
    stgs->acceleration_lookback = ACCELERATION_LOOKBACK;
    /* TODO add warm starting */
    stgs->warm_start = getIntFromListWithDefault(params, "warm_start", WARM_START);
    d->stgs = stgs;

    k->f = getIntFromListWithDefault(cone, "f", 0);
    k->l = getIntFromListWithDefault(cone, "l", 0);
    k->ep = getIntFromListWithDefault(cone, "ep", 0);
    k->ed = getIntFromListWithDefault(cone, "ed", 0);
    k->q = getIntVectorFromList(cone, "q", &(k->qsize));
    k->s = getIntVectorFromList(cone, "s", &(k->ssize));
    k->p = getFloatVectorFromList(cone, "p", &(k->psize));
    
    /* solve! */
    scs(d, k, sol, info);

    PROTECT(ret = NEW_LIST(4));
    PROTECT(retnames = NEW_CHARACTER(4));
    SET_NAMES(ret, retnames);
    
    SET_STRING_ELT(retnames, 0, mkChar("x"));
    SET_VECTOR_ELT(ret, 0, c_to_r_double_vec(sol->x, d->n));
    SET_STRING_ELT(retnames, 1, mkChar("y"));
    SET_VECTOR_ELT(ret, 1, c_to_r_double_vec(sol->y, d->m));
    SET_STRING_ELT(retnames, 2, mkChar("s"));
    SET_VECTOR_ELT(ret, 2, c_to_r_double_vec(sol->s, d->m));
    SET_STRING_ELT(retnames, 3, mkChar("info"));
    SET_VECTOR_ELT(ret, 3, info_to_r_list(info));
    
    /* free memory */
    scs_free(info);
    scs_free(d);
    scs_free(k);
    scs_free(stgs);
    scs_free(A);
    SCS(free_sol)(sol);
    UNPROTECT(2);
    return ret;
}


/**

void Rscs_finalize_data(ScsData *d) {
    scs_free(d);
}

void Rscs_finalize_cone(ScsCone *k) {
    scs_free(k);
}

void Rscs_finalize_settings(ScsSettings *s) {
    scs_free(s);
}

void Rscs_finalize_matrix(ScsMatrix *mat) {
    scs_free(mat);
}

void Rscs_finalize_info(ScsInfo *in) {
    scs_free(in);
}

void Rscs_finalize_solution(ScsSolution *solution) {
    free_sol(solution);
}

// Data
SEXP Rscs_new_data(void) {
    SEXP x, class_info;
    ScsData * d = scs_malloc(sizeof(ScsData));
    d->m = 0;
    d->n = 0;
    d->c = NULL;
    d->A = NULL;
    d->b = NULL;
    d->stgs = NULL;
    x = PROTECT(R_MakeExternalPtr(d, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(x, Rscs_finalize_data, TRUE);
    
    PROTECT(class_info = allocVector(STRSXP, 1));
    SET_STRING_ELT(class_info, 0, mkChar("scs_data"));
    setAttrib(x, install("class"), class_info);
    UNPROTECT(2);
    return x;
}

// ScsCone
SEXP Rscs_new_cone(void) {
    SEXP x, class_info;
    ScsCone * d = scs_malloc(sizeof(ScsCone));
    x = PROTECT(R_MakeExternalPtr(d, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(x, Rscs_finalize_cone, TRUE);
    
    PROTECT(class_info = allocVector(STRSXP, 1));
    SET_STRING_ELT(class_info, 0, mkChar("scs_cone"));
    setAttrib(x, install("class"), class_info);
    UNPROTECT(2);
    return x;
}

// ScsSettings
SEXP Rscs_new_settings(void) {
    SEXP x, class_info;
    ScsSettings * d = scs_malloc(sizeof(ScsSettings));
    x = PROTECT(R_MakeExternalPtr(d, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(x, Rscs_finalize_settings, TRUE);
    
    PROTECT(class_info = allocVector(STRSXP, 1));
    SET_STRING_ELT(class_info, 0, mkChar("scs_settings"));
    setAttrib(x, install("class"), class_info);
    UNPROTECT(2);
    return x;
}

// ScsMatrix
SEXP Rscs_new_matrix(void) {
    SEXP x, class_info;
    ScsMatrix * d = scs_malloc(sizeof(ScsMatrix));
    x = PROTECT(R_MakeExternalPtr(d, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(x, Rscs_finalize_data, TRUE);
    
    PROTECT(class_info = allocVector(STRSXP, 1));
    SET_STRING_ELT(class_info, 0, mkChar("scs_matrix"));
    setAttrib(x, install("class"), class_info);
    UNPROTECT(2);
    return x;
}

// ScsInfo
SEXP Rscs_new_info(void) {
    SEXP x, class_info;
    ScsInfo * d = scs_malloc(sizeof(ScsInfo));
    x = PROTECT(R_MakeExternalPtr(d, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(x, Rscs_finalize_data, TRUE);
    
    PROTECT(class_info = allocVector(STRSXP, 1));
    SET_STRING_ELT(class_info, 0, mkChar("scs_info"));
    setAttrib(x, install("class"), class_info);
    UNPROTECT(2);
    return x;
}

// ScsSolution
SEXP Rscs_new_solution(void) {
    SEXP x, class_info;
    ScsSolution * d = scs_malloc(sizeof(ScsSolution));
    x = PROTECT(R_MakeExternalPtr(d, R_NilValue, R_NilValue));
    R_RegisterCFinalizerEx(x, Rscs_finalize_solution, TRUE);
    
    PROTECT(class_info = allocVector(STRSXP, 1));
    SET_STRING_ELT(class_info, 0, mkChar("scs_data"));
    setAttrib(x, install("class"), class_info);
    UNPROTECT(2);
    return x;
}

SEXP Rscs_get_params(SEXP r_params) {
    SEXP r_list;
    ScsSettings *params = R_ExternalPtrAddr(r_params);
    
    PROTECT(r_list = allocVector(VECSXP, 8));
    
    SET_VECTOR_ELT(r_list, 0, c_to_r_integer(params->max_iters));
    SET_VECTOR_ELT(r_list, 1, c_to_r_integer(params->normalize));
    SET_VECTOR_ELT(r_list, 2, c_to_r_integer(params->verbose));
    
    SET_VECTOR_ELT(r_list, 3, c_to_r_double(params->cg_rate));
    SET_VECTOR_ELT(r_list, 4, c_to_r_double(params->scale));
    SET_VECTOR_ELT(r_list, 5, c_to_r_double(params->rho_x));
    SET_VECTOR_ELT(r_list, 6, c_to_r_double(params->alpha));
    SET_VECTOR_ELT(r_list, 7, c_to_r_double(params->eps));
    
    UNPROTECT(1);
    
    return r_list;
}

SEXP Rscs_set_params(SEXP r_params, SEXP max_iters, SEXP normalize, SEXP verbose,
                     SEXP cg_rate, SEXP scale, SEXP rho_x, SEXP alpha, SEXP eps) {
    ScsSettings *params = R_ExternalPtrAddr(r_params);
    
    params->max_iters  = INTEGER(max_iters)[0];
    params->normalize  = INTEGER(normalize)[0];
    params->verbose    = INTEGER(verbose)[0];
    params->cg_rate    = REAL(cg_rate)[0];
    params->scale      = REAL(scale)[0];
    params->rho_x      = REAL(rho_x)[0];
    params->alpha      = REAL(alpha)[0];
    params->eps        = REAL(eps)[0];
    
    return R_NilValue;
}

SEXP Rscs_data_get_matrix(SEXP r_dat) {
    ScsData *dat = R_ExternalPtrAddr(r_dat);
    if ( dat->A == NULL )
        return R_NilValue;
    
    ScsMatrix *A = dat->A;
    SEXP r_A;
    PROTECT(r_A = allocVector(VECSXP, 5));
    SET_VECTOR_ELT(r_A, 0, c_to_r_integer(A->m));
    SET_VECTOR_ELT(r_A, 1, c_to_r_integer(A->n));
    
    // TODO:
    
    return A;    
}

SEXP Rscs_data_set_matrix(SEXP r_dat, SEXP A_nrow, SEXP A_ncol, 
                          SEXP A_x, SEXP A_i, SEXP A_p) {
    ScsData *dat = R_ExternalPtrAddr(r_dat);
    if ( dat->A != NULL )
        free(A);
        
    ScsMatrix * A = scs_malloc(sizeof(ScsMatrix));
    A->m = INTEGER(A_nrow)[0];
    A->n = INTEGER(A_ncol)[0];
    A->x = REAL(A_x);
    A->i = INTEGER(A_i);
    A->p = INTEGER(A_p);
    
    return R_NilValue;
}

// returns null since anhow the values in the pointers changed
SEXP Rscs_solve(SEXP r_dat, SEXP r_cones, SEXP r_solution, SEXP r_info) {
    ScsData      *dat = R_ExternalPtrAddr(r_dat);
    ScsCone    *cones = R_ExternalPtrAddr(r_cones);
    ScsSolution  *solution = R_ExternalPtrAddr(r_solution);
    ScsSettings *info = R_ExternalPtrAddr(r_info);
      
    scs(dat, cones, solution, info);
    
    return R_NilValue;
}

*/

