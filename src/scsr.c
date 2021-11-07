#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include "scs/include/glbopts.h"
#include "scs/include/scs.h"
#include "scs/include/util.h"
#include "scs/linsys/scs_matrix.h"

SEXP getListElement(SEXP list, const char *str) {
  SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
  for (int i = 0; i < length(list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
      elmt = VECTOR_ELT(list, i);
      break;
    }
  }
  return elmt;
}

SEXP floatVec2R(scs_int n, scs_float *in) {
  SEXP ret;
  scs_int i;
  scs_float *vec;

  PROTECT(ret = allocVector(REALSXP, n));
  vec = REAL(ret);
  for (i = 0; i < n; i++) {
    vec[i] = in[i];
  }
  UNPROTECT(1);
  return ret;
}

const char *getStrFromListWithDefault(SEXP list, const char *str,
				      const char *def) {
  SEXP val = getListElement(list, str);
  if (val == R_NilValue) {
    return def;
  }
  val = coerceVector(val, STRSXP);
  return CHAR(STRING_ELT(val, 0));
}

scs_float getFloatFromListWithDefault(SEXP list, const char *str,
                                      scs_float def) {
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

scs_float *getFloatVectorFromList(SEXP list, const char *str, scs_int *len) {
  SEXP vec = getListElement(list, str);
  *len = length(vec);
  vec = coerceVector(vec, REALSXP);
  return REAL(vec);
}

scs_int *getIntVectorFromList(SEXP list, const char *str, scs_int *len) {
  SEXP vec = getListElement(list, str);
  *len = length(vec);
  vec = coerceVector(vec, INTSXP);
  return INTEGER(vec);
}

SEXP populateInfoR(ScsInfo *info) {
  scs_int num_protected = 0;
  SEXP infor, info_names, iter_r, status_r, statusVal_r, scaleUpdates_r, pobj_r, dobj_r,
    resPri_r, resDual_r, gap_r, resInfeas_r, resUnbddA_r, resUnbddP_r, 
    setupTime_r, solveTime_r, scale_r, compSlack_r, rejectedAccelSteps_r,
    acceptedAccelSteps_r, linsysTime_r, coneTime_r, accelTime_r;

  PROTECT(infor = NEW_LIST(21));
  num_protected++;

  PROTECT(info_names = NEW_CHARACTER(21));
  SET_NAMES(infor, info_names);
  UNPROTECT(1);
  
  PROTECT(iter_r = allocVector(INTSXP, 1));
  INTEGER(iter_r)[0] = info->iter;
  SET_STRING_ELT(info_names, 0, mkChar("iter"));
  SET_VECTOR_ELT(infor, 0, iter_r);
  UNPROTECT(1);

  
  PROTECT(status_r = NEW_CHARACTER(1));
  SET_STRING_ELT(status_r, 0, mkChar(info->status));
  SET_STRING_ELT(info_names, 1, mkChar("status"));
  SET_VECTOR_ELT(infor, 1, status_r);
  UNPROTECT(1);
  
  PROTECT(statusVal_r = allocVector(INTSXP, 1));
  INTEGER(statusVal_r)[0] = info->status_val;
  SET_STRING_ELT(info_names, 2, mkChar("status_val"));
  SET_VECTOR_ELT(infor, 2, statusVal_r);
  UNPROTECT(1);

  PROTECT(scaleUpdates_r = allocVector(INTSXP, 1));
  INTEGER(scaleUpdates_r)[0] = info->scale_updates;
  SET_STRING_ELT(info_names, 3, mkChar("scale_updates"));
  SET_VECTOR_ELT(infor, 3, scaleUpdates_r);
  UNPROTECT(1);
  
  PROTECT(pobj_r = allocVector(REALSXP, 1));
  REAL(pobj_r)[0] = info->pobj;
  SET_STRING_ELT(info_names, 4, mkChar("pobj"));
  SET_VECTOR_ELT(infor, 4, pobj_r);
  UNPROTECT(1);
  
  PROTECT(dobj_r = allocVector(REALSXP, 1));
  REAL(dobj_r)[0] = info->dobj;
  SET_STRING_ELT(info_names, 5, mkChar("dobj"));
  SET_VECTOR_ELT(infor, 5, dobj_r);
  UNPROTECT(1);
  
  PROTECT(resPri_r = allocVector(REALSXP, 1));
  REAL(resPri_r)[0] = info->res_pri;
  SET_STRING_ELT(info_names, 6, mkChar("res_pri"));
  SET_VECTOR_ELT(infor, 6, resPri_r);
  UNPROTECT(1);
  
  PROTECT(resDual_r = allocVector(REALSXP, 1));
  REAL(resDual_r)[0] = info->res_dual;
  SET_STRING_ELT(info_names, 7, mkChar("res_dual"));
  SET_VECTOR_ELT(infor, 7, resDual_r);
  UNPROTECT(1);

  PROTECT(gap_r = allocVector(REALSXP, 1));
  REAL(gap_r)[0] = info->gap;
  SET_STRING_ELT(info_names, 8, mkChar("gap"));
  SET_VECTOR_ELT(infor, 8, gap_r);
  UNPROTECT(1);
  
  PROTECT(resInfeas_r = allocVector(REALSXP, 1));
  REAL(resInfeas_r)[0] = info->res_infeas;
  SET_STRING_ELT(info_names, 9, mkChar("res_infeas"));
  SET_VECTOR_ELT(infor, 9, resInfeas_r);
  UNPROTECT(1);
  
  PROTECT(resUnbddA_r = allocVector(REALSXP, 1));
  REAL(resUnbddA_r)[0] = info->res_unbdd_a;
  SET_STRING_ELT(info_names, 10, mkChar("res_unbdd_a"));
  SET_VECTOR_ELT(infor, 10, resUnbddA_r);
  UNPROTECT(1);

  PROTECT(resUnbddP_r = allocVector(REALSXP, 1));
  REAL(resUnbddP_r)[0] = info->res_unbdd_p;
  SET_STRING_ELT(info_names, 11, mkChar("res_unbdd_p"));
  SET_VECTOR_ELT(infor, 11, resUnbddP_r);
  UNPROTECT(1);
  
  PROTECT(setupTime_r = allocVector(REALSXP, 1));
  REAL(setupTime_r)[0] = info->setup_time;
  SET_STRING_ELT(info_names, 12, mkChar("setup_time"));
  SET_VECTOR_ELT(infor, 12, setupTime_r);
  UNPROTECT(1);
  
  PROTECT(solveTime_r = allocVector(REALSXP, 1));
  REAL(solveTime_r)[0] = info->solve_time;
  SET_STRING_ELT(info_names, 13, mkChar("solve_time"));
  SET_VECTOR_ELT(infor, 13, solveTime_r);
  UNPROTECT(1);

  PROTECT(scale_r = allocVector(REALSXP, 1));
  REAL(scale_r)[0] = info->scale;
  SET_STRING_ELT(info_names, 14, mkChar("scale"));
  SET_VECTOR_ELT(infor, 14, scale_r);
  UNPROTECT(1);

  PROTECT(compSlack_r = allocVector(REALSXP, 1));
  REAL(compSlack_r)[0] = info->comp_slack;
  SET_STRING_ELT(info_names, 15, mkChar("comp_slack"));
  SET_VECTOR_ELT(infor, 15, compSlack_r);
  UNPROTECT(1);
  
  PROTECT(rejectedAccelSteps_r = allocVector(INTSXP, 1));
  INTEGER(rejectedAccelSteps_r)[0] = info->rejected_accel_steps;
  SET_STRING_ELT(info_names, 16, mkChar("rejected_accel_steps"));
  SET_VECTOR_ELT(infor, 16, rejectedAccelSteps_r);
  UNPROTECT(1);

  PROTECT(acceptedAccelSteps_r = allocVector(INTSXP, 1));
  INTEGER(acceptedAccelSteps_r)[0] = info->accepted_accel_steps;
  SET_STRING_ELT(info_names, 17, mkChar("accepted_accel_steps"));
  SET_VECTOR_ELT(infor, 17, acceptedAccelSteps_r);
  UNPROTECT(1);

  PROTECT(linsysTime_r = allocVector(REALSXP, 1));
  REAL(linsysTime_r)[0] = info->lin_sys_time;
  SET_STRING_ELT(info_names, 18, mkChar("lin_sys_time"));
  SET_VECTOR_ELT(infor, 18, linsysTime_r);
  UNPROTECT(1);

  PROTECT(coneTime_r = allocVector(REALSXP, 1));
  REAL(coneTime_r)[0] = info->cone_time;
  SET_STRING_ELT(info_names, 19, mkChar("cone_time"));
  SET_VECTOR_ELT(infor, 19, coneTime_r);
  UNPROTECT(1);
  
  PROTECT(accelTime_r = allocVector(REALSXP, 1));
  REAL(accelTime_r)[0] = info->accel_time;
  SET_STRING_ELT(info_names, 20, mkChar("accel_time"));
  SET_VECTOR_ELT(infor, 20, accelTime_r);
  UNPROTECT(1);

  UNPROTECT(num_protected); /* Unprotect the first PROTECT in this fn */
  return infor;
}

SEXP scsr(SEXP data, SEXP cone, SEXP params) {
  scs_int len;
  SEXP ret, retnames, infor, xr, yr, sr;

  /* allocate memory */
  ScsData *d = scs_malloc(sizeof(ScsData));
  ScsCone *k = scs_malloc(sizeof(ScsCone));
  ScsSettings *stgs = scs_malloc(sizeof(ScsSettings));
  ScsMatrix *A = scs_malloc(sizeof(ScsMatrix));
  ScsMatrix *P = SCS_NULL;
  ScsInfo *info = scs_calloc(1, sizeof(ScsInfo));
  ScsSolution *sol = scs_calloc(1, sizeof(ScsSolution));
  scs_int m, n;
  scs_float *fptr;
  
  d->b = getFloatVectorFromList(data, "b", &len);
  d->c = getFloatVectorFromList(data, "c", &len);
  n = d->n = getIntFromListWithDefault(data, "n", 0);
  m = d->m = getIntFromListWithDefault(data, "m", 0);

  A->m = m;
  A->n = n;
  A->x = getFloatVectorFromList(data, "Ax", &len);
  A->i = getIntVectorFromList(data, "Ai", &len);
  A->p = getIntVectorFromList(data, "Ap", &len);
  d->A = A;

  SEXP tmp = getListElement(data, "Px");
  if (tmp != R_NilValue) {
    P = scs_malloc(sizeof(ScsMatrix));
    P->x = getFloatVectorFromList(data, "Px", &len);
    P->i = getIntVectorFromList(data, "Pi", &len);
    P->p = getIntVectorFromList(data, "Pp", &len);
    P->m = n;
    P->n = n;
  }
  d->P = P;    

  /* New parameters for scs 3.0; see scs/src/include/glbopts.h for details*/
  stgs->max_iters = getIntFromListWithDefault(params, "max_iters", MAX_ITERS);
  stgs->eps_rel = getFloatFromListWithDefault(params, "eps_rel", EPS_REL);
  stgs->eps_abs = getFloatFromListWithDefault(params, "eps_abs", EPS_ABS);
  stgs->eps_infeas = getFloatFromListWithDefault(params, "eps_infeas", EPS_INFEAS);
  stgs->alpha = getFloatFromListWithDefault(params, "alpha", ALPHA);
  stgs->rho_x = getFloatFromListWithDefault(params, "rho_x", RHO_X);
  stgs->scale = getFloatFromListWithDefault(params, "scale", SCALE);
  stgs->verbose = getIntFromListWithDefault(params, "verbose", VERBOSE);
  stgs->normalize = getIntFromListWithDefault(params, "normalize", NORMALIZE);
  stgs->acceleration_lookback = getIntFromListWithDefault(params, "acceleration_lookback",
							  ACCELERATION_LOOKBACK);
  stgs->acceleration_interval = getIntFromListWithDefault(params, "acceleration_interval",
							  ACCELERATION_INTERVAL);
  stgs->adaptive_scale = getIntFromListWithDefault(params, "adaptive_scale", ADAPTIVE_SCALE);
  /* Without this nulling out, things bomb! */
  /* stgs->write_data_filename = NULL; */
  /* stgs->log_csv_filename = NULL; */
  stgs->write_data_filename = getStrFromListWithDefault(params, "write_data_filename", (char *) NULL);
  stgs->log_csv_filename = getStrFromListWithDefault(params, "log_csv_filename", (char *) NULL);
  stgs->time_limit_secs = getFloatFromListWithDefault(params, "time_limit_secs", TIME_LIMIT_SECS);

  /* Warm start data consists of x, y, s from previous solution and need to be stuffed into */
  /* the solution struct */
  SEXP initial = getListElement(data, "initial");
  if (initial == R_NilValue) {
    stgs->warm_start = WARM_START;
  } else {
    stgs->warm_start = 1;
    sol->x = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
    fptr = getFloatVectorFromList(initial, "x", &len);
    memcpy(sol->x, fptr, n * sizeof(scs_float));
    sol->y = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
    fptr = getFloatVectorFromList(initial, "y", &len);
    memcpy(sol->y, fptr, m * sizeof(scs_float));
    sol->s = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
    fptr = getFloatVectorFromList(initial, "s", &len);
    memcpy(sol->s, fptr, m * sizeof(scs_float));
  }
  
  k->z = getIntFromListWithDefault(cone, "z", 0);
  k->l = getIntFromListWithDefault(cone, "l", 0);
  k->ep = getIntFromListWithDefault(cone, "ep", 0);
  k->ed = getIntFromListWithDefault(cone, "ed", 0);
  k->qsize = getIntFromListWithDefault(cone, "", 0);
  k->q = getIntVectorFromList(cone, "q", &(k->qsize));
  k->s = getIntVectorFromList(cone, "s", &(k->ssize));
  k->p = getFloatVectorFromList(cone, "p", &(k->psize));
  k->bsize = getIntFromListWithDefault(cone, "bsize", 0);
  if (k->bsize > 0) {
    k->bl = getFloatVectorFromList(cone, "bl", &len);
    k->bu = getFloatVectorFromList(cone, "bu", &len);
  }

  /* solve! */
  /* scs(d, k, sol, info); */
  /* exit_flag is stuffed in info anyway, so ignore... */
  scs(d, k, stgs, sol, info);

  PROTECT(infor = populateInfoR(info)); /* count = 1 */

  PROTECT(ret = NEW_LIST(4));  /* count = 2 */

  PROTECT(retnames = NEW_CHARACTER(4));
  SET_NAMES(ret, retnames);
  UNPROTECT(1);
  
  PROTECT(xr = floatVec2R(d->n, sol->x));
  SET_STRING_ELT(retnames, 0, mkChar("x"));
  SET_VECTOR_ELT(ret, 0, xr);
  UNPROTECT(1);
  
  PROTECT(yr = floatVec2R(d->m, sol->y));
  SET_STRING_ELT(retnames, 1, mkChar("y"));
  SET_VECTOR_ELT(ret, 1, yr);
  UNPROTECT(1);

  PROTECT(sr = floatVec2R(d->m, sol->s));
  SET_STRING_ELT(retnames, 2, mkChar("s"));
  SET_VECTOR_ELT(ret, 2, sr);
  UNPROTECT(1);

  SET_STRING_ELT(retnames, 3, mkChar("info"));
  SET_VECTOR_ELT(ret, 3, infor);
  UNPROTECT(1);
  
  /* free memory */
  scs_free(info);
  scs_free(d);
  scs_free(k);
  scs_free(stgs);
  scs_free(A);
  if (tmp != R_NilValue) scs_free(P);
  SCS(free_sol)(sol);
  UNPROTECT(1);
  return ret;
}
