## Regenerate src/scs_mods/ from the pristine bundled SCS sources (src/scs).
##
## src/Makevars copies the files in src/scs_mods over the upstream tree before
## libscsdir.a is built; the submodule itself is never committed in modified
## form. This script is the only place the modifications are defined. Run it
## from the package root after every submodule bump:
##
##     Rscript tools/make_scs_mods.R
##
## Every edit names the exact upstream text it expects and how many times. If
## upstream changed that text the script stops with a DRIFT message instead of
## writing a stale copy. What is changed and why is listed in
## src/scs_mods/README.md; keep the two in step.

scs_root <- "src/scs"
mods_dir <- "src/scs_mods"
if (!dir.exists(scs_root) || !file.exists("DESCRIPTION")) {
  stop("run from the package root; '", scs_root, "' not found in ", getwd())
}
dirty <- suppressWarnings(system2("git", c("-C", scs_root, "status", "--short"),
                                  stdout = TRUE, stderr = FALSE))
if (length(dirty)) {
  stop("the submodule has local modifications (an in-tree build copies scs_mods ",
       "over it); restore it first:\n  git -C src/scs checkout -- . && ",
       "git -C src/scs clean -fdq\n", paste(dirty, collapse = "\n"))
}

count_fixed <- function(text, pattern) {
  lengths(regmatches(text, gregexpr(pattern, text, fixed = TRUE)))
}

## Replace the fixed string `from` by `to` in `text`, asserting that `from`
## occurs exactly `expect` times.
edit_fixed <- function(text, from, to, expect, label, file) {
  n <- count_fixed(text, from)
  if (n != expect) {
    stop("DRIFT [", label, "]: expected ", expect, " occurrence(s) of\n  ", from,
         "\nin ", file, " but found ", n,
         ". Upstream changed; review the modification before regenerating.")
  }
  cat(sprintf("%-4s %-12s %d site(s)\n", label, file, expect))
  gsub(from, to, text, fixed = TRUE)
}

read_upstream <- function(file) {
  path <- file.path(scs_root, file)
  if (!file.exists(path)) stop("DRIFT: missing ", path)
  paste(readLines(path, warn = FALSE), collapse = "\n")
}
write_mod <- function(text, file) {
  dir.create(mods_dir, showWarnings = FALSE)
  writeLines(text, file.path(mods_dir, basename(file)))
}

## 1. scs.mk: build flags ----------------------------------------------------
mk <- read_upstream("scs.mk")
mk <- edit_fixed(
  mk,
  from = "OPT ?= -O3 -march=native -fno-math-errno",
  to = paste("# R package: no optimization flags of its own, R's CFLAGS govern",
             "# (-march=native is not portable).",
             "OPT ?=", sep = "\n"),
  expect = 1L, label = "M1a", file = "scs.mk"
)
mk <- edit_fixed(
  mk,
  from = paste("override CFLAGS += -g -Wall -Wwrite-strings -pedantic",
               "-funroll-loops -Wstrict-prototypes $(INCLUDE) $(OPT)",
               "-Werror=incompatible-pointer-types"),
  to = "override CFLAGS += -g -Wall -pedantic -funroll-loops $(INCLUDE) $(OPT)",
  expect = 1L, label = "M1b", file = "scs.mk"
)
write_mod(mk, "scs.mk")

## 2. cones.c: sprintf -> snprintf ------------------------------------------
## 12 sites in format_cone_header(): 8 compiled always, 4 under
## USE_SPECTRAL_CONES (changed too; harmless). `line` is `char line[128]`.
cones <- read_upstream("src/cones.c")
cones <- edit_fixed(
  cones,
  from = "sprintf(line, ",
  to = "snprintf(line, sizeof(line), ",
  expect = 12L, label = "M2", file = "cones.c"
)

## 3. cones.c: complex PSD cone, zheevr -> zheev -----------------------------
## R's LAPACK has no zheevr. Three sites, each matched exactly.
cones <- edit_fixed(
  cones,
  from = paste(
    'void BLASC(heevr)(const char *jobz, const char *range, const char *uplo,',
    '                  blas_int *n, SCS_BLAS_COMPLEX_TYPE *a, blas_int *lda,',
    '                  scs_float *vl, scs_float *vu, blas_int *il, blas_int *iu,',
    '                  scs_float *abstol, blas_int *m, scs_float *w,',
    '                  SCS_BLAS_COMPLEX_TYPE *z, blas_int *ldz, blas_int *isuppz,',
    '                  SCS_BLAS_COMPLEX_TYPE *cwork, blas_int *lcwork,',
    '                  scs_float *rwork, blas_int *lrwork, blas_int *iwork,',
    '                  blas_int *liwork, blas_int *info);',
    sep = "\n"),
  to = paste(
    '/* R package: zheev instead of zheevr (R\'s LAPACK has no zheevr) */',
    'void BLASC(heev)(const char *jobz, const char *uplo, blas_int *n,',
    '                 SCS_BLAS_COMPLEX_TYPE *a, blas_int *lda, scs_float *w,',
    '                 SCS_BLAS_COMPLEX_TYPE *work, blas_int *lwork,',
    '                 scs_float *rwork, blas_int *info);',
    sep = "\n"),
  expect = 1L, label = "M3a", file = "cones.c"
)
cones <- edit_fixed(
  cones,
  from = paste(
    '    BLASC(heevr)("V", "A", "L", &n_max_csd, SCS_BLAS_COMPLEX_CAST(c->cXs),',
    '                 &n_max_csd, &d_f, &d_f, &d_i, &d_i, &abstol, &m, c->e,',
    '                 SCS_BLAS_COMPLEX_CAST(c->cZ), &n_max_csd, c->isuppz,',
    '                 SCS_BLAS_COMPLEX_CAST(&lcwork_opt), &neg_one, &lrwork_opt,',
    '                 &neg_one, &liwork_opt_c, &neg_one, &info);',
    sep = "\n"),
  to = paste(
    '    /* R package: zheev workspace query; rwork is not touched by a query */',
    '    BLASC(heev)("V", "L", &n_max_csd, SCS_BLAS_COMPLEX_CAST(c->cXs),',
    '                &n_max_csd, c->e, SCS_BLAS_COMPLEX_CAST(&lcwork_opt),',
    '                &neg_one, &lrwork_opt, &info);',
    '    /* zheev needs rwork of length max(1, 3n - 2); it lives in c->work */',
    '    lrwork_opt = (scs_float)MAX(1, 3 * n_max_csd - 2);',
    sep = "\n"),
  expect = 1L, label = "M3b", file = "cones.c"
)
cones <- edit_fixed(
  cones,
  from = paste(
    '  BLASC(heevr)("V", "A", "L", &nb, SCS_BLAS_COMPLEX_CAST(c->cXs), &nb, &d_f,',
    '               &d_f, &d_i, &d_i, &abstol, &m, c->e,',
    '               SCS_BLAS_COMPLEX_CAST(c->cZ), &nb, c->isuppz,',
    '               SCS_BLAS_COMPLEX_CAST(c->cwork), &c->lcwork, c->work, &c->lwork,',
    '               c->iwork, &c->liwork, &info);',
    sep = "\n"),
  to = paste(
    '  /* R package: zheev (R\'s LAPACK has no zheevr). The eigenvectors come',
    '   * back in cXs; copy them to cZ where the reconstruction expects them. */',
    '  (void)d_f; (void)d_i; (void)abstol; (void)m;',
    '  BLASC(heev)("V", "L", &nb, SCS_BLAS_COMPLEX_CAST(c->cXs), &nb, c->e,',
    '              SCS_BLAS_COMPLEX_CAST(c->cwork), &c->lcwork, c->work, &info);',
    '  if (info == 0)',
    '    memcpy(c->cZ, c->cXs, (size_t)n * (size_t)n * sizeof(scs_complex_float));',
    sep = "\n"),
  expect = 1L, label = "M3c", file = "cones.c"
)
write_mod(cones, "src/cones.c")

## Guard: no bare printf in any file the static archive compiles -------------
## (output must go through scs_printf, which is Rprintf under R_LANG)
compiled <- c(
  Sys.glob(file.path(scs_root, "src", "*.c")),
  Sys.glob(file.path(scs_root, "include", "*.h")),
  Sys.glob(file.path(scs_root, "linsys", "*.c")),
  Sys.glob(file.path(scs_root, "linsys", "*.h")),
  Sys.glob(file.path(scs_root, "linsys", "cpu", "direct", "*.[ch]")),
  Sys.glob(file.path(scs_root, "linsys", "external", "amd", "*.[ch]")),
  Sys.glob(file.path(scs_root, "linsys", "external", "qdldl", "*.[ch]"))
)
if (length(compiled) < 30L) {
  stop("DRIFT [guard]: only ", length(compiled), " compiled source files under ", scs_root)
}
bare_printf <- function(path) {
  lines <- readLines(path, warn = FALSE)
  hits <- grepl("(^|[^A-Za-z_])printf\\(", lines) & !grepl("scs_printf", lines)
  if (any(hits)) paste0(path, ":", which(hits)) else character()
}
hits <- unlist(lapply(compiled, bare_printf), use.names = FALSE)
if (length(hits)) {
  stop("DRIFT [guard]: bare printf( in compiled sources:\n  ", paste(hits, collapse = "\n  "))
}
cat(sprintf("%-4s bare printf guard clean (%d files)\n", "G", length(compiled)))

## Show the result as a diff against upstream (review aid; not stored) -------
for (f in c("scs.mk", "src/cones.c")) {
  d <- suppressWarnings(system2("diff", c("-u", file.path(scs_root, f),
                                          file.path(mods_dir, basename(f))),
                                stdout = TRUE, stderr = FALSE))
  cat(sprintf("-- %s: %d added, %d removed lines vs upstream\n", basename(f),
              sum(startsWith(d, "+")) - 1L, sum(startsWith(d, "-")) - 1L))
}
cat("make_scs_mods.R: src/scs_mods regenerated\n")
