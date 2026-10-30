/* Run the C test problems bundled with SCS (src/scs/test/problems) against a
 * library built with the R package's flags and dump, for every solve, the
 * exact problem data, cone, settings, warm start and the solution/info that
 * SCS returned.  dumps_to_rds.R turns each dump into one RDS reference file
 * for inst/tinytest/problem_results.  Driven by regenerate.sh.
 *
 * Mechanism: the test headers call scs() or scs_init()/scs_solve(); those
 * names are macro-redirected to wrappers that dump and then call the real
 * functions.  The wrappers are defined BEFORE the macros so they call the
 * real library.  Built as an R shared object (R CMD SHLIB) and run INSIDE an
 * R process via .Call("dump_tests_run", outdir), so the library sees exactly
 * the environment the package sees (same process, same BLAS/LAPACK bindings,
 * same Rprintf).
 */
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>
#include "scs.h"
#include "glbopts.h"

#define printf Rprintf

int tests_run = 0;
static const char *current_test = "";
static const char *out_dir = ".";
static int case_no = 0;
static int cur_solve_no = 1;   /* 1 for scs(); k for the k-th scs_solve() on one workspace */
static FILE *manifest = NULL;
/* a copy of the solution gen_random_prob_data() planted in the current test,
 * if any: upstream asserts c'x and b'y against c'x* and b'y* (small_lp.h:37-38).
 * Copied because the test owns and may free its opt_sol before a later solve. */
static scs_float *planted_x = NULL, *planted_y = NULL;
static scs_int planted_n = 0, planted_m = 0;
static void planted_clear(void) {
  free(planted_x); free(planted_y);
  planted_x = planted_y = NULL; planted_n = planted_m = 0;
}

static void wvec_f(FILE *f, const char *tag, const scs_float *v, scs_int n) {
  scs_int i;
  fprintf(f, "%s", tag);
  if (v) for (i = 0; i < n; i++) fprintf(f, " %a", v[i]);
  fprintf(f, "\n");
}
static void wvec_i(FILE *f, const char *tag, const scs_int *v, scs_int n) {
  scs_int i;
  fprintf(f, "%s", tag);
  if (v) for (i = 0; i < n; i++) fprintf(f, " %ld", (long)v[i]);
  fprintf(f, "\n");
}

static FILE *open_dump(const ScsData *d, const ScsCone *k, const ScsSettings *s,
                       scs_int warm_start, const ScsSolution *sol) {
  char path[1024];
  FILE *f;
  scs_int m = d ? d->m : -1, n = d ? d->n : -1;
  int valid = d && d->A && d->A->p && m > 0 && n > 0 && d->b && d->c;
  case_no++;
  snprintf(path, sizeof path, "%s/%s__%d.txt", out_dir, current_test, case_no);
  f = fopen(path, "w");
  if (!f) Rf_error("cannot open %s", path);
  fprintf(f, "test %s\ncase %d\nvalid %d\nm %ld\nn %ld\n", current_test, case_no, valid, (long)m, (long)n);
  /* A's own dims (upstream validation corrupts them; a dgCMatrix cannot) and
   * the solve count on this workspace (re-solves keep AA/diag_r state: not
   * reproducible by the one-shot R API) */
  fprintf(f, "A_m %ld\nA_n %ld\nsolve_no %d\n", (long)(d && d->A ? d->A->m : -1), (long)(d && d->A ? d->A->n : -1), cur_solve_no);
  if (valid) {
    scs_int nnz = d->A->p[n];
    fprintf(f, "A %ld\n", (long)nnz);
    wvec_i(f, "Ai", d->A->i, nnz); wvec_i(f, "Ap", d->A->p, n + 1); wvec_f(f, "Ax", d->A->x, nnz);
    if (d->P && d->P->p) {
      scs_int pnz = d->P->p[n];
      fprintf(f, "P %ld\n", (long)pnz);
      wvec_i(f, "Pi", d->P->i, pnz); wvec_i(f, "Pp", d->P->p, n + 1); wvec_f(f, "Px", d->P->x, pnz);
    } else {
      fprintf(f, "P -1\n");
    }
    wvec_f(f, "b", d->b, m); wvec_f(f, "c", d->c, n);
    if (planted_x && planted_y && planted_n == n && planted_m == m) {
      scs_float cx = 0, by = 0;
      scs_int i;
      for (i = 0; i < n; i++) cx += d->c[i] * planted_x[i];
      for (i = 0; i < m; i++) by += d->b[i] * planted_y[i];
      fprintf(f, "planted %a %a\n", cx, by);
    }
  }
  if (k) {
    fprintf(f, "cone_scalars %ld %ld %ld %ld %ld %ld %ld %ld %ld\n", (long)k->z, (long)k->l, (long)k->bsize,
            (long)k->ep, (long)k->ed, (long)k->qsize, (long)k->ssize, (long)k->cssize, (long)k->psize);
    wvec_f(f, "bl", k->bl, k->bsize > 0 ? k->bsize - 1 : 0); wvec_f(f, "bu", k->bu, k->bsize > 0 ? k->bsize - 1 : 0);
    wvec_i(f, "q", k->q, k->qsize > 0 ? k->qsize : 0); wvec_i(f, "s", k->s, k->ssize > 0 ? k->ssize : 0);
    wvec_i(f, "cs", k->cs, k->cssize > 0 ? k->cssize : 0); wvec_f(f, "p", k->p, k->psize > 0 ? k->psize : 0);
  } else {
    fprintf(f, "cone_scalars NA\n");
  }
  if (s) {
    fprintf(f, "settings %ld %a %ld %ld %a %ld %a %a %a %a %a %ld %ld %ld %ld %ld %a %a\n",
            (long)s->normalize, s->scale, (long)s->adaptive_scale, (long)s->adaptive_diag_scale, s->rho_x,
            (long)s->max_iters, s->eps_abs, s->eps_rel, s->eps_infeas, s->alpha, s->time_limit_secs,
            (long)s->verbose, (long)s->warm_start, (long)s->acceleration_lookback, (long)s->acceleration_interval,
            (long)s->acceleration_type_1, s->acceleration_regularization, s->acceleration_relaxation);
    fprintf(f, "write_data_filename %s\n", s->write_data_filename ? s->write_data_filename : "-");
  } else {
    fprintf(f, "settings NA\n");
  }
  fprintf(f, "warm_start %ld\n", (long)warm_start);
  if (warm_start && valid && sol) {
    wvec_f(f, "init_x", sol->x, n); wvec_f(f, "init_y", sol->y, m); wvec_f(f, "init_s", sol->s, m);
  }
  return f;
}

static void close_dump(FILE *f, scs_int exitflag, const ScsData *d, const ScsSolution *sol, const ScsInfo *info) {
  int valid = d && d->m > 0 && d->n > 0;
  fprintf(f, "exitflag %ld\n", (long)exitflag);
  if (info) {
    fprintf(f, "info %ld %ld %ld %a %a %a %a %a %a %a %a %a %a %ld %ld\n",
            (long)info->iter, (long)info->status_val, (long)info->scale_updates, info->pobj, info->dobj,
            info->res_pri, info->res_dual, info->gap, info->res_infeas, info->res_unbdd_a, info->res_unbdd_p,
            info->scale, info->comp_slack, (long)info->rejected_accel_steps, (long)info->accepted_accel_steps);
    fprintf(f, "status %s\nlin_sys_solver %s\n", info->status, info->lin_sys_solver);
    fprintf(f, "aa_stats %ld %ld %ld %ld %ld %ld %ld %ld %a %a\n", (long)info->aa_stats.iter,
            (long)info->aa_stats.n_accept, (long)info->aa_stats.n_reject_lapack, (long)info->aa_stats.n_reject_rank0,
            (long)info->aa_stats.n_reject_nonfinite, (long)info->aa_stats.n_reject_weight_cap,
            (long)info->aa_stats.n_safeguard_reject, (long)info->aa_stats.last_rank, info->aa_stats.last_aa_norm,
            info->aa_stats.last_regularization);
  }
  if (valid && sol && sol->x && sol->y && sol->s) {
    wvec_f(f, "sol_x", sol->x, d->n); wvec_f(f, "sol_y", sol->y, d->m); wvec_f(f, "sol_s", sol->s, d->m);
  }
  fprintf(f, "end\n");
  fclose(f);
}

/* ---- wrappers (call the real library) ---- */
static scs_int dump_scs(const ScsData *d, const ScsCone *k, const ScsSettings *s, ScsSolution *sol, ScsInfo *info) {
  FILE *f;
  cur_solve_no = 1;
  f = open_dump(d, k, s, s ? s->warm_start : 0, sol);
  scs_int r = scs(d, k, s, sol, info);
  close_dump(f, r, d, sol, info);
  return r;
}

static const ScsData *init_d = NULL;
static const ScsCone *init_k = NULL;
static const ScsSettings *init_s = NULL;
static int init_solves = 0;
static ScsWork *dump_scs_init(const ScsData *d, const ScsCone *k, const ScsSettings *s) {
  init_d = d; init_k = k; init_s = s; init_solves = 0;
  return scs_init(d, k, s);
}
static scs_int dump_scs_solve(ScsWork *w, ScsSolution *sol, ScsInfo *info, scs_int warm_start) {
  FILE *f;
  cur_solve_no = ++init_solves;
  f = open_dump(init_d, init_k, init_s, warm_start, sol);
  scs_int r = scs_solve(w, sol, info, warm_start);
  close_dump(f, r, init_d, sol, info);
  return r;
}

/* problem_utils.h is include-guarded: pulling it in here defines the real
 * generator once, so the macro below redirects only the tests' calls */
#include "problem_utils.h"
static void dump_gen_random_prob_data(scs_int nnz, scs_int col_nnz, ScsData *d, ScsCone *k,
                                      ScsSolution *opt_sol, scs_int seed) {
  gen_random_prob_data(nnz, col_nnz, d, k, opt_sol, seed);
  planted_clear();
  planted_n = d->n; planted_m = d->m;
  planted_x = (scs_float *)malloc(planted_n * sizeof(scs_float));
  planted_y = (scs_float *)malloc(planted_m * sizeof(scs_float));
  if (!planted_x || !planted_y) Rf_error("out of memory copying the planted solution");
  memcpy(planted_x, opt_sol->x, planted_n * sizeof(scs_float));
  memcpy(planted_y, opt_sol->y, planted_m * sizeof(scs_float));
}

#define gen_random_prob_data(nnz, col_nnz, d, k, opt_sol, seed) \
  dump_gen_random_prob_data(nnz, col_nnz, d, k, opt_sol, seed)
#define scs(d, k, s, sol, info) dump_scs(d, k, s, sol, info)
#define scs_init(d, k, s) dump_scs_init(d, k, s)
#define scs_solve(w, sol, info, ws) dump_scs_solve(w, sol, info, ws)

/* upstream test headers (spectral-cone, read/write, scs_update and internal
 * function tests are left out: not reachable through the R API) */
#include "problems/complex_PSD.h"
#include "problems/degenerate.h"
#include "problems/dense_qp.h"
#include "problems/hs21_tiny_qp.h"
#include "problems/infeasible_lp.h"
#include "problems/infeasible_socp.h"
#include "problems/infeasible_tiny_qp.h"
#include "problems/max_ent.h"
#include "problems/mpc_bug.h"
#include "problems/qafiro_tiny_qp.h"
#include "problems/random_prob.h"
#include "problems/rob_gauss_cov_est.h"
#include "problems/sd_and_complex_sd.h"
#include "problems/small_lp.h"
#include "problems/small_qp.h"
#include "problems/test_box_cone.h"
#include "problems/test_dual_exp_cone.h"
#include "problems/test_exp_cone.h"
#include "problems/test_inaccurate.h"
#include "problems/test_mixed_cones.h"
#include "problems/test_power_cone.h"
#include "problems/test_psd_metric.h"
#include "problems/test_psd_n1.h"
#include "problems/test_soc_sizes.h"
#include "problems/test_solver_options.h"
#include "problems/test_validation.h"
#include "problems/test_zero_cone.h"
#include "problems/unbounded_lp.h"
#include "problems/unbounded_socp.h"
#include "problems/unbounded_tiny_qp.h"

#define RUN(t)                                                        \
  do {                                                                \
    const char *msg;                                                  \
    current_test = #t; case_no = 0; planted_clear();                  \
    printf("******** %s\n", #t);                                      \
    msg = t();                                                        \
    tests_run++;                                                      \
    fprintf(manifest, "%s\t%d\t%s\n", #t, case_no, msg ? msg : "PASS"); \
    if (msg) { printf("UPSTREAM ASSERTION FAILED in %s: %s\n", #t, msg); failures++; } \
  } while (0)

SEXP dump_tests_run(SEXP outdir) {
  int failures = 0;
  char path[1024];
  SEXP res;
  tests_run = 0;
  out_dir = CHAR(STRING_ELT(outdir, 0));
  snprintf(path, sizeof path, "%s/MANIFEST.tsv", out_dir);
  manifest = fopen(path, "w");
  if (!manifest) Rf_error("cannot open %s", path);
  RUN(test_validation);
  RUN(degenerate);
  RUN(dense_qp);
  RUN(small_lp);
  RUN(small_qp);
  RUN(rob_gauss_cov_est);
  RUN(complex_PSD);
  RUN(sd_and_complex_sd);
  RUN(test_psd_metric);
  RUN(hs21_tiny_qp);
  RUN(qafiro_tiny_qp);
  RUN(infeasible_tiny_qp);
  RUN(infeasible_lp);
  RUN(infeasible_socp);
  RUN(unbounded_tiny_qp);
  RUN(unbounded_lp);
  RUN(unbounded_socp);
  RUN(random_prob);
  RUN(max_ent);
  RUN(mpc_bug);
  RUN(test_exp_cone);
  RUN(test_dual_exp_cone);
  RUN(test_power_cone);
  RUN(test_power_cone_p09);
  RUN(test_dual_power_cone);
  RUN(test_multi_power);
  RUN(test_power_cone_infeasible);
  RUN(test_soc_size1);
  RUN(test_soc_size2);
  RUN(test_soc_size3);
  RUN(test_soc_size5);
  RUN(test_multi_soc);
  RUN(test_zero_cone);
  RUN(test_box_cone_lp);
  RUN(test_psd_n1);
  RUN(test_solved_inaccurate);
  RUN(test_infeasible_inaccurate);
  RUN(test_unbounded_inaccurate);
  RUN(test_max_iters_1);
  RUN(test_adaptive_scale);
  RUN(test_no_acceleration);
  RUN(test_type2_acceleration);
  RUN(test_aa_relaxation_sweep);
  RUN(test_aa_regularization_sweep);
  RUN(test_negative_lookback_rejected);
  RUN(test_invalid_aa_relaxation_rejected);
  RUN(test_pinned_negative_regularization_solves);
  RUN(test_invalid_aa_regularization_rejected);
  RUN(test_normalize_off);
  RUN(test_warm_start);
  RUN(test_mixed_cones);
  fclose(manifest);
  manifest = NULL;
  planted_clear();
  printf("\n%d upstream tests run, %d upstream assertion failures, %s\n", tests_run, failures, out_dir);
  res = PROTECT(Rf_allocVector(INTSXP, 2));
  INTEGER(res)[0] = tests_run;
  INTEGER(res)[1] = failures;
  UNPROTECT(1);
  return res;
}
