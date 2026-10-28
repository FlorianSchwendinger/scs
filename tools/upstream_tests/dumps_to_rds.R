## Convert the dumps written by dump_tests (upstream SCS 3.3.1 test problems,
## solved by the oracle library) into the RDS layout used by
## inst/tinytest/test_problems.R:
##   list(data = list(A, P, b, obj), cone, settings, initial?,
##        upstream = list(test, case, scs, exitflag, planted?))
## A is a dgCMatrix and P a dsCMatrix (upper), exactly as the old files.
## Run: Rscript dumps_to_rds.R <dumps dir> <out dir>
args <- commandArgs(trailingOnly = TRUE)
dumps <- args[1]; out <- args[2]
dir.create(out, showWarnings = FALSE, recursive = TRUE)
suppressMessages(library(Matrix))

parse_dump <- function(path) {
  lines <- readLines(path)
  key <- sub(" .*$", "", lines)
  val <- sub("^[^ ]*( |$)", "", lines)
  get <- function(k) val[match(k, key)]
  ## Doubles are dumped as C99 hex floats (%a): exact, unlike R's parsing of
  ## 17-digit decimals, which is off by one ulp for a fifth of all values.
  hex2num <- function(tok) {
    tok <- tolower(tok)
    tok[tok %in% c("nan", "-nan")] <- "NaN"; tok[tok == "na"] <- "NA"
    tok[tok == "inf"] <- "Inf"; tok[tok == "-inf"] <- "-Inf"
    out <- suppressWarnings(as.numeric(tok))   ## the check below is the real guard
    bad <- is.na(out) & !(tok %in% c("NaN", "NA"))
    if (any(bad)) stop("unparsed numeric tokens: ", paste(head(tok[bad]), collapse = " "))
    out
  }
  num <- function(k) { v <- get(k); if (is.na(v) || !nzchar(v)) numeric(0) else hex2num(strsplit(v, " ", fixed = TRUE)[[1]]) }
  int <- function(k) as.integer(num(k))
  list(test = get("test"), case = as.integer(get("case")), valid = as.integer(get("valid")),
       m = as.integer(get("m")), n = as.integer(get("n")),
       A_m = as.integer(get("A_m")), A_n = as.integer(get("A_n")), solve_no = as.integer(get("solve_no")),
       A = list(i = int("Ai"), p = int("Ap"), x = num("Ax")),
       P = if (identical(get("P"), "-1") || is.na(get("P"))) NULL else list(i = int("Pi"), p = int("Pp"), x = num("Px")),
       b = num("b"), c = num("c"), planted = num("planted"),
       write_data_filename = get("write_data_filename"),
       cone_scalars = num("cone_scalars"), bl = num("bl"), bu = num("bu"), q = int("q"), s = int("s"), cs = int("cs"), p = num("p"),
       settings = num("settings"), warm_start = as.integer(get("warm_start")),
       init = if (!is.na(match("init_x", key))) list(x = num("init_x"), y = num("init_y"), s = num("init_s")) else NULL,
       exitflag = as.integer(get("exitflag")), info = num("info"), status = get("status"), lin_sys_solver = get("lin_sys_solver"),
       aa_stats = num("aa_stats"), sol_x = num("sol_x"), sol_y = num("sol_y"), sol_s = num("sol_s"))
}

settings_names <- c("normalize", "scale", "adaptive_scale", "adaptive_diag_scale", "rho_x", "max_iters", "eps_abs", "eps_rel",
                    "eps_infeas", "alpha", "time_limit_secs", "verbose", "warm_start", "acceleration_lookback",
                    "acceleration_interval", "acceleration_type_1", "acceleration_regularization", "acceleration_relaxation")
settings_int <- c("normalize", "adaptive_scale", "adaptive_diag_scale", "max_iters", "verbose", "warm_start",
                  "acceleration_lookback", "acceleration_interval", "acceleration_type_1")

to_rds <- function(d) {
  n <- d$n; m <- d$m
  A <- new("dgCMatrix", i = d$A$i, p = d$A$p, x = d$A$x, Dim = c(m, n))
  P <- if (is.null(d$P)) NULL else new("dsCMatrix", i = d$P$i, p = d$P$p, x = d$P$x, Dim = c(n, n), uplo = "U")
  cs <- d$cone_scalars
  cone <- list(z = as.integer(cs[1]), l = as.integer(cs[2]), bsize = as.integer(cs[3]), bl = d$bl, bu = d$bu,
               q = d$q, s = d$s, cs = d$cs, ep = as.integer(cs[4]), ed = as.integer(cs[5]), p = d$p)
  settings <- setNames(as.list(d$settings), settings_names)
  for (k in settings_int) settings[[k]] <- as.integer(settings[[k]])
  settings$verbose <- 0L   ## quiet tests; printing does not affect the iterates
  ## the solution SCS returned is not stored: test_problems.R checks the
  ## replayed solution with upstream's assertions, not against a saved one
  res <- list(data = list(A = A, P = P, b = d$b, obj = d$c), cone = cone, settings = settings,
              upstream = list(test = d$test, case = d$case, scs = "3.3.1", exitflag = d$exitflag))
  ## c'x* and b'y* of the solution gen_random_prob_data() planted: upstream's
  ## objective check for generated problems (small_lp.h:37-38)
  if (length(d$planted) == 2L) res$upstream$planted <- c(cx = d$planted[1], by = d$planted[2])
  if (d$warm_start == 1L && !is.null(d$init)) res$initial <- d$init
  res
}

files <- list.files(dumps, pattern = "__[0-9]+\\.txt$", full.names = TRUE)
parsed <- lapply(files, parse_dump)
tests <- vapply(parsed, `[[`, "", "test")
cases <- vapply(parsed, `[[`, 1L, "case")
ncase <- table(tests)
skipped <- character(0); written <- character(0)
for (i in order(tests, cases)) {
  d <- parsed[[i]]
  nm <- if (ncase[[d$test]] > 1L) paste0(d$test, d$case) else d$test
  reproducible <- d$valid == 1L && length(d$sol_x) > 0L && length(d$settings) == 18L && length(d$cone_scalars) == 9L
  if (!reproducible) { skipped <- c(skipped, sprintf("%s (valid=%d, settings=%d, cone=%d, exitflag=%d)", nm, d$valid, length(d$settings), length(d$cone_scalars), d$exitflag)); next }
  ## upstream corrupts A's own dims in one validation case; a dgCMatrix has one Dim
  if (d$A_m != d$m || d$A_n != d$n) { skipped <- c(skipped, sprintf("%s (A dims %dx%d disagree with data %dx%d: not expressible)", nm, d$A_m, d$A_n, d$m, d$n)); next }
  ## a k-th scs_solve() on one workspace keeps the AA memory and diag_r of the
  ## previous solve (aa_reset only runs on a scale update); the R API is one-shot
  if (d$solve_no > 1L) { skipped <- c(skipped, sprintf("%s (re-solve #%d on a reused workspace: not a one-shot computation)", nm, d$solve_no)); next }
  ## a call that only writes the problem to a data file (rob_gauss_cov_est.h:195-198)
  ## asserts nothing upstream, and the R API has no data-file writer
  if (!is.na(d$write_data_filename) && d$write_data_filename != "-") { skipped <- c(skipped, sprintf("%s (writes data file '%s'; upstream asserts nothing on this call)", nm, d$write_data_filename)); next }
  ## upstream validation cases with malformed matrices cannot exist as Matrix objects
  obj <- tryCatch(to_rds(d), error = function(e) { skipped <<- c(skipped, sprintf("%s (not constructible: %s)", nm, conditionMessage(e))); NULL })
  if (is.null(obj)) next
  saveRDS(obj, file.path(out, paste0(nm, ".RDS")), version = 2)
  written <- c(written, nm)
}
cat(length(written), "RDS files written to", out, "\n")
cat("skipped (no reproducible problem data):\n  ", paste(skipped, collapse = "\n   "), "\n")
