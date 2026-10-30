/*
 * R binding for the bundled SCS library (src/scs, upstream tag 3.3.1).
 *
 * Entry points (registered in init.c):
 *   scsr(data, cone, params)  solve one problem, return list(x, y, s, info)
 *   scsr_version()            SCS_VERSION string of the bundled library
 *
 * Conventions:
 *   - Every SEXP obtained from coerceVector() is PROTECTed for the whole call
 *     (counter `np`), because the C structs keep raw pointers into those
 *     vectors until scs() returns.
 *   - Settings start from scs_set_default_settings() and only the entries
 *     present in `params` override them, so a new upstream field can never be
 *     left uninitialized.
 *   - Cone and info fields follow the order of the ScsCone / ScsInfo structs
 *     in src/scs/include/scs.h.
 *   - Input shape errors are caught in R (scs()) before .Call; the checks
 *     here are a second line of defense and use error(), which longjmps.
 */

#include <string.h>
#include <R.h>
#include <Rinternals.h>

#include "scs/include/glbopts.h"
#include "scs/include/scs.h"
#include "scs/include/util.h" /* SCS(free_sol) */

/* ---------------------------------------------------------------- helpers */

static SEXP list_get(SEXP list, const char *name) {
  SEXP names;
  R_xlen_t i;
  if (list == R_NilValue || !isNewList(list)) return R_NilValue;
  names = getAttrib(list, R_NamesSymbol);
  if (names == R_NilValue) return R_NilValue;
  for (i = 0; i < xlength(list); i++) {
    if (strcmp(CHAR(STRING_ELT(names, i)), name) == 0) {
      return VECTOR_ELT(list, i);
    }
  }
  return R_NilValue;
}

/* Coerce and PROTECT; the caller unprotects *np at the end of scsr(). */
static SEXP coerce_keep(SEXP x, SEXPTYPE type, int *np) {
  SEXP y = PROTECT(coerceVector(x, type));
  (*np)++;
  return y;
}

static scs_int get_int(SEXP list, const char *name, scs_int def, int *np) {
  SEXP v = list_get(list, name);
  if (v == R_NilValue || xlength(v) == 0) return def;
  v = coerce_keep(v, INTSXP, np);
  if (INTEGER(v)[0] == NA_INTEGER) error("'%s' must not be NA", name);
  return (scs_int)INTEGER(v)[0];
}

static scs_float get_real(SEXP list, const char *name, scs_float def,
                          int *np) {
  SEXP v = list_get(list, name);
  if (v == R_NilValue || xlength(v) == 0) return def;
  v = coerce_keep(v, REALSXP, np);
  return (scs_float)REAL(v)[0];
}

/* NULL when absent, empty, or NA. */
static const char *get_string(SEXP list, const char *name, int *np) {
  SEXP v = list_get(list, name);
  if (v == R_NilValue || xlength(v) == 0) return SCS_NULL;
  v = coerce_keep(v, STRSXP, np);
  if (STRING_ELT(v, 0) == NA_STRING) return SCS_NULL;
  return CHAR(STRING_ELT(v, 0));
}

/* Vectors: NULL with *len = 0 when absent or empty. */
static scs_float *get_real_vec(SEXP list, const char *name, scs_int *len,
                               int *np) {
  SEXP v = list_get(list, name);
  *len = 0;
  if (v == R_NilValue || xlength(v) == 0) return SCS_NULL;
  v = coerce_keep(v, REALSXP, np);
  *len = (scs_int)xlength(v);
  return REAL(v);
}

static scs_int *get_int_vec(SEXP list, const char *name, scs_int *len,
                            int *np) {
  SEXP v = list_get(list, name);
  *len = 0;
  if (v == R_NilValue || xlength(v) == 0) return SCS_NULL;
  v = coerce_keep(v, INTSXP, np);
  *len = (scs_int)xlength(v);
  return INTEGER(v);
}

/* Matrix arrays: never NULL. SCS 3.3.1 rejects NULL x/i/p even when a matrix
   has no stored entries (e.g. an all-zero P), so an empty array points at a
   dummy that is never read (nnz = 0). */
static scs_float empty_real[1];
static scs_int empty_int[1];

static scs_float *get_mat_real(SEXP list, const char *name, scs_int *len,
                               int *np) {
  scs_float *x = get_real_vec(list, name, len, np);
  return x ? x : empty_real;
}

static scs_int *get_mat_int(SEXP list, const char *name, scs_int *len,
                            int *np) {
  scs_int *x = get_int_vec(list, name, len, np);
  return x ? x : empty_int;
}

static void copy_vec(const char *what, SEXP list, const char *name,
                     scs_float *dst, scs_int n, int *np) {
  scs_int len;
  scs_float *src = get_real_vec(list, name, &len, np);
  if (len != n) {
    error("'%s$%s' must have length %ld, not %ld", what, name, (long)n,
          (long)len);
  }
  memcpy(dst, src, (size_t)n * sizeof(scs_float));
}

/* Fresh REALSXP from a C array; NA-filled when the array is missing. */
static SEXP real_from_c(const scs_float *x, scs_int n) {
  SEXP r = PROTECT(allocVector(REALSXP, n));
  scs_int i;
  for (i = 0; i < n; i++) REAL(r)[i] = x ? x[i] : NA_REAL;
  UNPROTECT(1);
  return r;
}

/* ------------------------------------------------------------- info list */

typedef enum { F_INT, F_REAL, F_STR, F_SEXP } field_kind;

typedef struct {
  const char *name;
  field_kind kind;
  const void *ptr; /* scs_int*, scs_float*, const char*, or SEXP* */
} field_t;

static SEXP build_list(const field_t *f, int n) {
  SEXP out = PROTECT(allocVector(VECSXP, n));
  SEXP names = PROTECT(allocVector(STRSXP, n));
  int i;
  for (i = 0; i < n; i++) {
    SEXP v;
    switch (f[i].kind) {
    case F_INT:
      v = ScalarInteger((int)*(const scs_int *)f[i].ptr);
      break;
    case F_REAL:
      v = ScalarReal((double)*(const scs_float *)f[i].ptr);
      break;
    case F_STR:
      v = mkString((const char *)f[i].ptr);
      break;
    default: /* F_SEXP: an already protected SEXP */
      v = *(SEXP *)f[i].ptr;
      break;
    }
    SET_VECTOR_ELT(out, i, v);
    SET_STRING_ELT(names, i, mkChar(f[i].name));
  }
  setAttrib(out, R_NamesSymbol, names);
  UNPROTECT(2);
  return out;
}

static SEXP aa_stats_to_R(const AaStats *a) {
  const field_t f[] = {
      {"iter", F_INT, &a->iter},
      {"n_accept", F_INT, &a->n_accept},
      {"n_reject_lapack", F_INT, &a->n_reject_lapack},
      {"n_reject_rank0", F_INT, &a->n_reject_rank0},
      {"n_reject_nonfinite", F_INT, &a->n_reject_nonfinite},
      {"n_reject_weight_cap", F_INT, &a->n_reject_weight_cap},
      {"n_safeguard_reject", F_INT, &a->n_safeguard_reject},
      {"last_rank", F_INT, &a->last_rank},
      {"last_aa_norm", F_REAL, &a->last_aa_norm},
      {"last_regularization", F_REAL, &a->last_regularization},
  };
  return build_list(f, (int)(sizeof(f) / sizeof(f[0])));
}

static SEXP info_to_R(const ScsInfo *info) {
  SEXP aa = PROTECT(aa_stats_to_R(&info->aa_stats));
  const field_t f[] = {
      {"iter", F_INT, &info->iter},
      {"status", F_STR, info->status},
      {"lin_sys_solver", F_STR, info->lin_sys_solver},
      {"status_val", F_INT, &info->status_val},
      {"scale_updates", F_INT, &info->scale_updates},
      {"pobj", F_REAL, &info->pobj},
      {"dobj", F_REAL, &info->dobj},
      {"res_pri", F_REAL, &info->res_pri},
      {"res_dual", F_REAL, &info->res_dual},
      {"gap", F_REAL, &info->gap},
      {"res_infeas", F_REAL, &info->res_infeas},
      {"res_unbdd_a", F_REAL, &info->res_unbdd_a},
      {"res_unbdd_p", F_REAL, &info->res_unbdd_p},
      {"setup_time", F_REAL, &info->setup_time},
      {"solve_time", F_REAL, &info->solve_time},
      {"scale", F_REAL, &info->scale},
      {"comp_slack", F_REAL, &info->comp_slack},
      {"rejected_accel_steps", F_INT, &info->rejected_accel_steps},
      {"accepted_accel_steps", F_INT, &info->accepted_accel_steps},
      {"aa_stats", F_SEXP, &aa},
      {"lin_sys_time", F_REAL, &info->lin_sys_time},
      {"cone_time", F_REAL, &info->cone_time},
      {"accel_time", F_REAL, &info->accel_time},
  };
  SEXP out = build_list(f, (int)(sizeof(f) / sizeof(f[0])));
  UNPROTECT(1);
  return out;
}

/* --------------------------------------------------------------- settings */

static void read_settings(ScsSettings *s, SEXP params, int *np) {
  scs_set_default_settings(s);
  /* Order follows ScsSettings in src/scs/include/scs.h. */
  s->normalize = get_int(params, "normalize", s->normalize, np);
  s->scale = get_real(params, "scale", s->scale, np);
  s->adaptive_scale = get_int(params, "adaptive_scale", s->adaptive_scale, np);
  s->adaptive_diag_scale =
      get_int(params, "adaptive_diag_scale", s->adaptive_diag_scale, np);
  s->rho_x = get_real(params, "rho_x", s->rho_x, np);
  s->max_iters = get_int(params, "max_iters", s->max_iters, np);
  s->eps_abs = get_real(params, "eps_abs", s->eps_abs, np);
  s->eps_rel = get_real(params, "eps_rel", s->eps_rel, np);
  s->eps_infeas = get_real(params, "eps_infeas", s->eps_infeas, np);
  s->alpha = get_real(params, "alpha", s->alpha, np);
  s->time_limit_secs =
      get_real(params, "time_limit_secs", s->time_limit_secs, np);
  s->verbose = get_int(params, "verbose", s->verbose, np);
  /* warm_start is set from `initial` in scsr(), not from params. */
  s->acceleration_lookback =
      get_int(params, "acceleration_lookback", s->acceleration_lookback, np);
  s->acceleration_interval =
      get_int(params, "acceleration_interval", s->acceleration_interval, np);
  s->acceleration_type_1 =
      get_int(params, "acceleration_type_1", s->acceleration_type_1, np);
  s->acceleration_regularization = get_real(
      params, "acceleration_regularization", s->acceleration_regularization,
      np);
  s->acceleration_relaxation = get_real(
      params, "acceleration_relaxation", s->acceleration_relaxation, np);
  s->write_data_filename = get_string(params, "write_data_filename", np);
  s->log_csv_filename = get_string(params, "log_csv_filename", np);
}

/* ------------------------------------------------------------------ cones */

static void read_cone(ScsCone *k, SEXP cone, int *np) {
  scs_int len_bl, len_bu;
  /* k comes from scs_calloc: every pointer NULL, every size 0. */
  k->z = get_int(cone, "z", 0, np);
  k->l = get_int(cone, "l", 0, np);
  k->bsize = get_int(cone, "bsize", 0, np);
  if (k->bsize > 0) {
    k->bl = get_real_vec(cone, "bl", &len_bl, np);
    k->bu = get_real_vec(cone, "bu", &len_bu, np);
    if (len_bl != k->bsize - 1 || len_bu != k->bsize - 1) {
      error("box cone: 'bl' and 'bu' must have length bsize - 1 = %ld",
            (long)(k->bsize - 1));
    }
  }
  k->q = get_int_vec(cone, "q", &k->qsize, np);
  k->s = get_int_vec(cone, "s", &k->ssize, np);
  k->cs = get_int_vec(cone, "cs", &k->cssize, np);
  k->ep = get_int(cone, "ep", 0, np);
  k->ed = get_int(cone, "ed", 0, np);
  k->p = get_real_vec(cone, "p", &k->psize, np);
}

/* ------------------------------------------------------------------- data */

static void read_data(ScsData *d, ScsMatrix *A, ScsMatrix **P, SEXP data,
                      int *np) {
  scs_int len, nnz;
  d->n = get_int(data, "n", 0, np);
  d->m = get_int(data, "m", 0, np);
  if (d->m <= 0 || d->n <= 0) {
    error("the problem must have at least one row and one column");
  }
  d->b = get_real_vec(data, "b", &len, np);
  if (len != d->m) error("'b' must have length nrow(A) = %ld", (long)d->m);
  d->c = get_real_vec(data, "c", &len, np);
  if (len != d->n) error("'obj' must have length ncol(A) = %ld", (long)d->n);

  A->m = d->m;
  A->n = d->n;
  A->x = get_mat_real(data, "Ax", &nnz, np);
  A->i = get_mat_int(data, "Ai", &len, np);
  if (len != nnz) error("'Ai' and 'Ax' must have the same length");
  A->p = get_mat_int(data, "Ap", &len, np);
  if (len != d->n + 1) error("'Ap' must have length ncol(A) + 1");
  d->A = A;

  *P = SCS_NULL;
  if (list_get(data, "Px") != R_NilValue) {
    ScsMatrix *Pm = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
    Pm->m = d->n;
    Pm->n = d->n;
    Pm->x = get_mat_real(data, "Px", &nnz, np);
    Pm->i = get_mat_int(data, "Pi", &len, np);
    if (len != nnz) error("'Pi' and 'Px' must have the same length");
    Pm->p = get_mat_int(data, "Pp", &len, np);
    if (len != d->n + 1) error("'Pp' must have length ncol(A) + 1");
    *P = Pm;
  }
  d->P = *P;
}

/* ----------------------------------------------------------- entry points */

SEXP scsr(SEXP data, SEXP cone, SEXP params) {
  int np = 0; /* number of PROTECTs held for the duration of the call */
  ScsData *d = (ScsData *)scs_calloc(1, sizeof(ScsData));
  ScsCone *k = (ScsCone *)scs_calloc(1, sizeof(ScsCone));
  ScsSettings *stgs = (ScsSettings *)scs_calloc(1, sizeof(ScsSettings));
  ScsMatrix *A = (ScsMatrix *)scs_calloc(1, sizeof(ScsMatrix));
  ScsMatrix *P = SCS_NULL;
  ScsInfo *info = (ScsInfo *)scs_calloc(1, sizeof(ScsInfo));
  ScsSolution *sol = (ScsSolution *)scs_calloc(1, sizeof(ScsSolution));
  SEXP initial, x, y, s, infor, ret;

  read_data(d, A, &P, data, &np);
  read_settings(stgs, params, &np);
  read_cone(k, cone, &np);

  initial = list_get(data, "initial");
  if (initial == R_NilValue) {
    stgs->warm_start = 0; /* SCS allocates x, y, s */
  } else {
    stgs->warm_start = 1;
    sol->x = (scs_float *)scs_calloc(d->n, sizeof(scs_float));
    sol->y = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
    sol->s = (scs_float *)scs_calloc(d->m, sizeof(scs_float));
    copy_vec("initial", initial, "x", sol->x, d->n, &np);
    copy_vec("initial", initial, "y", sol->y, d->m, &np);
    copy_vec("initial", initial, "s", sol->s, d->m, &np);
  }

  /* The exit flag is also info->status_val; on failure SCS fills x, y, s
     with NaN itself (populate_on_failure). */
  scs(d, k, stgs, sol, info);

  x = PROTECT(real_from_c(sol->x, d->n));
  y = PROTECT(real_from_c(sol->y, d->m));
  s = PROTECT(real_from_c(sol->s, d->m));
  infor = PROTECT(info_to_R(info));
  {
    const field_t f[] = {
        {"x", F_SEXP, &x},
        {"y", F_SEXP, &y},
        {"s", F_SEXP, &s},
        {"info", F_SEXP, &infor},
    };
    ret = PROTECT(build_list(f, 4));
  }

  SCS(free_sol)(sol);
  scs_free(info);
  scs_free(d);
  scs_free(k);
  scs_free(stgs);
  scs_free(A);
  if (P) scs_free(P);

  UNPROTECT(5 + np);
  return ret;
}

SEXP scsr_version(void) { return mkString(scs_version()); }
