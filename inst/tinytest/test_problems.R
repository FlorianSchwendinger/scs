## "Saved Problems"
##
## Every file in problem_results/ is one solve from the SCS 3.3.2 C test suite
## (src/scs/test/problems/*.h; see problem_results/README.md). Each one is
## replayed through scs() and checked with the assertions upstream makes for
## that solve: the exit flag, the objective against upstream's known optimum,
## and verify_solution_correct() (test/problem_utils.h), all with upstream's
## tolerances. These are properties of the returned solution, so they hold on
## every platform.

## ------------------------------------------------------------ linear algebra
## Matrix objects are read through their slots only.
csc_mult <- function(M, v) {             # M %*% v
    out <- numeric(M@Dim[1L])
    if (length(M@x)) {
        j <- rep.int(seq_len(M@Dim[2L]), diff(M@p))
        r <- rowsum(M@x * v[j], M@i + 1L)
        out[as.integer(rownames(r))] <- r[, 1L]
    }
    out
}
csc_tmult <- function(M, v) {            # t(M) %*% v
    out <- numeric(M@Dim[2L])
    if (length(M@x)) {
        j <- rep.int(seq_len(M@Dim[2L]), diff(M@p))
        r <- rowsum(M@x * v[M@i + 1L], j)
        out[as.integer(rownames(r))] <- r[, 1L]
    }
    out
}
sym_upper_mult <- function(P, v) {       # P %*% v, P symmetric with its upper triangle stored
    j <- rep.int(seq_len(P@Dim[2L]), diff(P@p))
    d <- (P@i + 1L) == j
    diag_part <- numeric(length(v))
    diag_part[j[d]] <- P@x[d] * v[j[d]]
    csc_mult(P, v) + csc_tmult(P, v) - diag_part
}
norm_inf <- function(v) if (length(v)) max(abs(v)) else 0     # NORM is SCS(norm_inf), glbopts.h:227

## ---------------------------------------------------------- cone projections
## Projections onto the cone K for the zero, nonnegative, box and second-order
## cones (cones.c). PSD, complex PSD, exponential and power blocks are not
## projected: their rows are left out of the cone-distance checks, every other
## check below still applies to them.
proj_box <- function(tx, bl, bu) {      # cones.c:1825-1888, unscaled (r_box = NULL), t warm start 1
    if (length(tx) == 1L) return(max(tx, 0))
    x <- tx[-1L]; t <- 1
    above <- function(t) { h <- t * bu; !is.na(h) & x > h }
    below <- function(t) { l <- t * bl; !is.na(l) & x < l }
    for (iter in seq_len(25L)) {         # BOX_CONE_MAX_ITERS
        t_prev <- t
        a <- above(t); b <- below(t) & !a
        gt <- (t - tx[1L]) + sum((t * bu[a] - x[a]) * bu[a]) + sum((t * bl[b] - x[b]) * bl[b])
        ht <- 1 + sum(bu[a]^2) + sum(bl[b]^2)
        t <- max(t - gt / max(ht, 1e-8), 0)
        if (abs(gt / max(ht, 1e-6)) < 1e-12 * max(t, 1) || abs(t - t_prev) < 1e-11 * max(t, 1)) break
    }
    a <- above(t); b <- below(t) & !a
    x[a] <- t * bu[a]; x[b] <- t * bl[b]
    c(t, x)
}
proj_soc <- function(x) {               # cones.c proj_soc
    v1 <- x[1L]; s <- sqrt(sum(x[-1L]^2))
    if (s <= v1) return(x)
    if (s <= -v1) return(numeric(length(x)))
    alpha <- (v1 + s) / 2
    c(alpha, alpha * x[-1L] / s)
}
## Pi_K(v) on the projected rows, and which rows those are
proj_cone <- function(v, cone) {
    out <- v; checked <- logical(length(v)); pos <- 0L
    block <- function(len, f) {
        if (len > 0L) {
            r <- pos + seq_len(len)
            if (!is.null(f)) { out[r] <<- f(v[r]); checked[r] <<- TRUE }
            pos <<- pos + len
        }
    }
    num <- function(x) if (is.null(x)) 0L else as.integer(x)
    block(num(cone$z), function(v) numeric(length(v)))
    block(num(cone$l), function(v) pmax(v, 0))
    block(num(cone$bsize), function(v) proj_box(v, cone$bl, cone$bu))
    for (q in cone$q) block(as.integer(q), proj_soc)
    for (k in cone$s) block(as.integer(k * (k + 1) / 2), NULL)
    for (k in cone$cs) block(as.integer(k * k), NULL)
    block(3L * num(cone$ep), NULL)
    block(3L * num(cone$ed), NULL)
    block(3L * length(cone$p), NULL)
    list(v = out, checked = checked)
}
## ||s - Pi_K(s)|| = ||Pi_K*(-s)||, and ||Pi_K*(y) - y||, where by Moreau
## Pi_K*(v) = v + Pi_K(-v) (problem_utils.h:168-190)
pri_cone_dist <- function(s, cone) { p <- proj_cone(s, cone); norm_inf((p$v - s)[p$checked]) }
dual_cone_dist <- function(y, cone) { p <- proj_cone(-y, cone); norm_inf(p$v[p$checked]) }

## -------------------------------------------- verify_solution_correct(), R port
## test/problem_utils.h:192-352, double precision (VERIFY_TOL takes its first
## argument). mu_assert_less(msg, a, b) fails when a > b.
SOLVED <- 1L; UNBOUNDED <- -1L; INFEASIBLE <- -2L
verify_solution_correct <- function(d, sol, status, fname) {
    less <- function(msg, a, b) expect_false(a > b, info = paste(fname, msg, a, ">", b))
    A <- d$data$A; P <- d$data$P; b <- d$data$b; c <- d$data$obj
    stgs <- d$settings; info <- sol$info
    x <- sol$x; y <- sol$y; s <- sol$s
    ax <- csc_mult(A, x)
    primal <- ax + s
    res_unbdd_a <- norm_inf(primal)
    res_pri <- norm_inf(primal - b)
    if (!is.null(P)) {
        px <- sym_upper_mult(P, x); xt_p_x <- sum(px * x); res_unbdd_p <- norm_inf(px)
    } else {
        px <- numeric(length(x)); xt_p_x <- 0; res_unbdd_p <- 0
    }
    aty <- csc_tmult(A, y)
    res_infeas <- norm_inf(aty)
    res_dual <- norm_inf(aty + px + c)
    sdist <- if (status %in% c(SOLVED, UNBOUNDED)) pri_cone_dist(s, d$cone) else NaN
    ydist <- if (status %in% c(SOLVED, INFEASIBLE)) dual_cone_dist(y, d$cone) else NaN
    sty <- sum(y * s); bty <- sum(y * b); ctx <- sum(x * c)
    gap <- abs(xt_p_x + ctx + bty)
    pobj <- xt_p_x / 2 + ctx
    dobj <- -xt_p_x / 2 - bty
    grl <- max(abs(xt_p_x), abs(ctx), abs(bty))
    prl <- max(norm_inf(b), norm_inf(s), norm_inf(ax))
    drl <- max(norm_inf(c), norm_inf(px), norm_inf(aty))
    if (status == SOLVED) {
        less("Primal residual ERROR", abs(res_pri - info$res_pri), 1e-10 * (1 + abs(pobj)))
        less("Dual residual ERROR", abs(res_dual - info$res_dual), 1e-10 * (1 + abs(pobj)))
        less("Gap ERROR", abs(gap - info$gap), 1e-7 * (1 + abs(gap)))
        less("Primal obj ERROR", abs(pobj - info$pobj), 1e-9 * (1 + abs(pobj)))
        less("Dual obj ERROR", abs(dobj - info$dobj), 1e-9 * (1 + abs(dobj)))
        less("Complementary slackness ERROR", abs(sty), 5e-8 * max(norm_inf(s), norm_inf(y)))
        less("s cone dist ERROR", abs(sdist), 1e-5)
        less("y cone dist ERROR", abs(ydist), 1e-5)
        less("Primal feas ERROR", res_pri, stgs$eps_abs + stgs$eps_rel * prl)
        less("Dual feas ERROR", res_dual, stgs$eps_abs + stgs$eps_rel * drl)
        less("Gap feas ERROR", gap, stgs$eps_abs + stgs$eps_rel * grl)
    } else if (status == INFEASIBLE) {
        less("Infeas ERROR", abs(res_infeas - info$res_infeas), 1e-8)
        less("bty ERROR", abs(bty + 1), 1e-12)
        less("y cone dist ERROR", abs(ydist), 1e-5)
        less("Infeas invalid ERROR", res_infeas, stgs$eps_infeas)
    } else if (status == UNBOUNDED) {
        less("Unbdd_a ERROR", abs(res_unbdd_a - info$res_unbdd_a), 1e-8)
        less("Unbdd_p ERROR", abs(res_unbdd_p - info$res_unbdd_p), 1e-8)
        less("ctx ERROR", abs(ctx + 1), 1e-12)
        less("s cone dist ERROR", abs(sdist), 1e-5)
        less("Unbounded P invalid ERROR", res_unbdd_p, stgs$eps_infeas)
        less("Unbounded A invalid ERROR", res_unbdd_a, stgs$eps_infeas)
    } else {
        expect_true(FALSE, info = paste(fname, "INVALID STATUS"))
    }
}

## ---------------------------------------------------- upstream's assertions
## flag: the exit flag upstream asserts. opt/form: upstream's objective check,
## |perr| < obj_tol && |derr| < obj_tol (SOL_TOL = 1e-4 unless the test sets
## its own), with the test's own formula:
##   "info"    perr = info$pobj - opt,            derr = info$dobj - opt
##   "cxby"    perr = c'x - opt,                  derr = -b'y - opt
##   "quad"    perr = x'Px/2 + c'x - opt,         derr = -x'Px/2 - b'y - opt
##   "planted" perr = c'x - c'x*,                 derr = -b'y + b'y*
##             (x*, y*: the solution gen_random_prob_data() planted; recorded
##              in upstream$planted)
## verify: upstream calls verify_solution_correct() after this solve.
## extra: any further upstream assertion; `res` holds the solves of this run.
u <- function(flag, opt = NULL, form = "info", verify = FALSE, extra = NULL, obj_tol = 1e-4)
    list(flag = flag, opt = opt, form = form, verify = verify, extra = extra, obj_tol = obj_tol)
same_obj <- function(ref, tol) function(sol, d, res) {      # degenerate.h:92-96
    expect_true(abs(sol$info$pobj - res[[ref]]$info$pobj) < tol &&
                abs(sol$info$dobj - res[[ref]]$info$dobj) < tol)
}
iter_le_half <- function(ref) function(sol, d, res)        # warm start, relative to the cold solve
    expect_true(sol$info$iter <= res[[ref]]$info$iter %/% 2L)
iter_le <- function(cap) function(sol, d, res) expect_true(sol$info$iter <= cap)
iter_is_max <- function(sol, d, res) expect_identical(sol$info$iter, as.integer(d$settings$max_iters))
iter_is_1 <- function(sol, d, res) expect_identical(sol$info$iter, 1L)
## the validation cases: SCS rejects the input and returns SCS_FAILED
validation <- u(-4L)

expected <- list(
    complex_PSD = u(1L, -5.228930, "cxby", TRUE),                       # complex_PSD.h:62-74
    degenerate1 = u(1L, verify = TRUE),                                 # degenerate.h:73-79
    degenerate2 = u(1L, verify = TRUE, extra = same_obj("degenerate1", 1e-8)),   # :87-96
    degenerate3 = u(1L, verify = TRUE, extra = same_obj("degenerate1", 1e-8)),   # :103-112
    degenerate4 = u(1L, verify = TRUE, extra = same_obj("degenerate1", 1e-8)),   # :120-129
    dense_qp = u(1L, form = "planted", verify = TRUE),                  # dense_qp.h:45-57
    hs21_tiny_qp1 = u(1L, 0.04000000000000625, "info", TRUE),           # hs21_tiny_qp.h:74-87
    ## warm start: upstream re-tests call 1's perr/derr here (not recomputed)
    hs21_tiny_qp2 = u(1L, extra = iter_le_half("hs21_tiny_qp1")),      # :91-100
    hs21_tiny_qp3 = u(1L, 0.04000000000000625, "info"),                 # :104-123
    infeasible_lp = u(-2L, verify = TRUE),                              # infeasible_lp.h:54-58
    infeasible_socp = u(-2L, verify = TRUE),                            # infeasible_socp.h:63-67
    infeasible_tiny_qp = u(-2L, verify = TRUE),                         # infeasible_tiny_qp.h:87-92
    issue_140 = u(1L, 2.875113, "quad", TRUE),                          # issue_140.h:16-19, test_prob_from_data_file.h:42-65
    issue_220 = u(1L, 3825 / 4096, "quad", TRUE, obj_tol = 5e-3),       # issue_220.h:34-37, test_prob_from_data_file.h:42-65
    issue_499 = u(-2L, verify = TRUE),                                  # issue_499.h:67-71
    max_ent = u(1L, -6.067087663361563, "quad", TRUE),                  # max_ent.h:6, test_prob_from_data_file.h:42-65
    mpc_bug1 = u(1L, -0.473957794500, "quad", TRUE),                    # mpc_bug.h:7-18
    mpc_bug2 = u(1L, -0.029336830816, "quad", TRUE),
    mpc_bug3 = u(1L, -0.002215217478, "quad", TRUE),
    qafiro_tiny_qp1 = u(1L, -1.5907818, "info", TRUE),                  # qafiro_tiny_qp.h:129-140
    qafiro_tiny_qp2 = u(1L, -1.5907818, "info", extra = iter_le(100L)), # :142-160, warm start
    qafiro_tiny_qp3 = u(1L, -1.5907818, "info"),                        # :167-187
    qafiro_tiny_qp4 = u(1L, -1.5907818, "info"),                        # :189-211
    random_prob = u(1L, 5.751458006385587, "quad", TRUE),               # random_prob.h:6
    rob_gauss_cov_est1 = u(1L, -4.8912, "cxby", TRUE),                  # rob_gauss_cov_est.h:180-194
    ## warm start: upstream re-tests call 1's perr/derr here (not recomputed)
    rob_gauss_cov_est2 = u(1L, extra = iter_le_half("rob_gauss_cov_est1")),  # :200-208
    rob_gauss_cov_est4 = u(1L, -4.8912, "info", TRUE),                  # :228-243
    rob_gauss_cov_est5 = u(1L, -4.8912, "info"),                        # :246-267
    sd_and_complex_sd = u(1L, -135.009872, "info", TRUE),               # sd_and_complex_sd.h:129-143
    small_lp = u(1L, form = "planted", verify = TRUE),                  # small_lp.h:35-47
    small_qp1 = u(1L, verify = TRUE),                                   # small_qp.h
    test_aa_regularization_sweep1 = u(1L, verify = TRUE),               # test_solver_options.h:252-259
    test_aa_regularization_sweep2 = u(1L, verify = TRUE),
    test_aa_regularization_sweep3 = u(1L, verify = TRUE),
    test_aa_relaxation_sweep1 = u(1L, verify = TRUE),                   # test_solver_options.h:218-225
    test_aa_relaxation_sweep2 = u(1L, verify = TRUE),
    test_aa_relaxation_sweep3 = u(1L, verify = TRUE),
    test_adaptive_scale = u(1L, verify = TRUE, extra = function(sol, d, res) {   # test_solver_options.h:97-105
        expect_true(abs(sum(d$data$obj * sol$x) - d$upstream$planted[["cx"]]) < 1e-3)
        expect_true(sol$info$scale_updates > 0L)
    }),
    test_box_cone_lp = u(1L, -3.0, "info", TRUE),                       # test_box_cone.h:76-85
    test_dual_exp_cone = u(1L, 0.13533528323661270, "info", TRUE),      # test_dual_exp_cone.h:80-89
    test_dual_power_cone = u(1L, -2.0, "info", TRUE),                   # test_power_cone.h:243-252
    test_infeasible_inaccurate = u(-7L, extra = iter_is_max),           # test_inaccurate.h:111-116
    test_invalid_aa_regularization_rejected1 = validation,              # test_solver_options.h:332-339
    test_invalid_aa_regularization_rejected2 = validation,
    test_invalid_aa_regularization_rejected3 = validation,
    test_invalid_aa_relaxation_rejected = validation,                   # test_solver_options.h:278-284
    test_max_iters_1 = u(2L, extra = iter_is_1),                        # test_inaccurate.h:195-200
    test_mixed_cones = u(1L, 4.0, "info", TRUE),                        # test_mixed_cones.h:97-106
    test_multi_power = u(1L, -3.0, "info", TRUE),                       # test_power_cone.h:401-410
    test_multi_soc = u(1L, 3.41421356237309504, "info", TRUE),          # test_soc_sizes.h:381-390
    ## this package rejects a negative lookback in R, before SCS sees it
    test_negative_lookback_rejected = c(u(-4L), r_error = TRUE),        # test_solver_options.h:187-194
    test_no_acceleration = u(1L, verify = TRUE, extra = function(sol, d, res) {  # test_solver_options.h:131-139
        expect_identical(sol$info$accepted_accel_steps, 0L)
        expect_true(sol$info$aa_stats$iter == 0L && sol$info$aa_stats$n_accept == 0L &&
                    sol$info$aa_stats$n_safeguard_reject == 0L)
    }),
    ## this build's linear solver (QDLDL) checks the inertia of the KKT matrix,
    ## so upstream asserts SCS_FAILED; an error raised in R would not count
    test_non_psd_p = u(-4L),                                            # test_non_psd_p.h:62-67
    test_normalize_off1 = u(1L, verify = TRUE),                         # test_solver_options.h:361-368
    test_normalize_off2 = u(1L, verify = TRUE, extra = function(sol, d, res)     # :371-381
        expect_true(abs(sol$info$pobj - res[["test_normalize_off1"]]$info$pobj) < 1e-4)),
    test_pinned_negative_regularization_solves = u(1L),                 # test_solver_options.h:304-311
    test_power_cone_infeasible = u(-2L, verify = TRUE),                 # test_power_cone.h:316-320
    test_power_cone_p09 = u(1L, -1.0, "info", TRUE),                    # test_power_cone.h:158-167
    test_power_cone_small_scale = u(1L, verify = TRUE, extra = function(sol, d, res) {  # test_power_cone.h:484-497
        viol <- abs(sol$s[3L]) - max(sol$s[1L], 0)^0.25 * max(sol$s[2L], 0)^0.75
        expect_true(abs(sol$x[1L] - 0.1) < 1e-6 && viol < 1e-6)
    }),
    test_power_cone1 = u(1L, -2.0, "info", TRUE),                       # test_power_cone.h:70-79
    test_power_cone2 = u(1L, extra = function(sol, d, res) {            # :81-86, warm start
        expect_true(sol$info$iter <= 25L)
        expect_true(abs(sol$info$pobj - (-2.0)) < 1e-4)
    }),
    test_psd_n1 = u(1L, 1.0, "info", TRUE),                             # test_psd_n1.h:68-77
    test_psd_scale_updates = u(1L, verify = TRUE, extra = function(sol, d, res) {   # test_psd_scale_updates.h:161-173
        expect_false(abs(sol$info$pobj - 235.52015) > 1e-2)
        expect_true(sol$info$scale_updates <= 10L || sol$info$scale_updates * 250L <= sol$info$iter)
    }),
    test_soc_size1 = u(1L, 1.5, "info", TRUE),                          # test_soc_sizes.h:63-72
    test_soc_size2 = u(1L, 1.5, "info", TRUE),                          # :138-147
    test_soc_size3 = u(1L, 2.2360679774997896, "info", TRUE),           # :214-223
    test_soc_size5 = u(1L, 2.0, "info", TRUE),                          # :295-304
    test_solved_inaccurate = u(2L, extra = iter_is_max),                # test_inaccurate.h:66-72
    test_type2_acceleration = u(1L, verify = TRUE),                     # test_solver_options.h:159-168
    test_unbounded_inaccurate = u(-6L, extra = iter_is_max),            # test_inaccurate.h:155-160
    test_validation1 = validation, test_validation2 = validation,       # test_validation.h:51-271
    test_validation3 = validation, test_validation4 = validation,
    test_validation5 = validation, test_validation6 = validation,
    test_validation7 = validation, test_validation8 = validation,
    test_validation9 = validation, test_validation10 = validation,
    test_validation11 = validation, test_validation12 = validation,
    test_validation15 = validation, test_validation22 = validation,
    test_warm_start1 = u(1L, verify = TRUE),                            # test_solver_options.h:536-544
    test_zero_cone = u(1L, 3.0, "info", TRUE),                          # test_zero_cone.h:65-74
    unbounded_lp = u(-1L, verify = TRUE),                               # unbounded_lp.h:54-58
    unbounded_socp = u(-1L, verify = TRUE),                             # unbounded_socp.h:59-63
    unbounded_tiny_qp = u(-1L, verify = TRUE)                           # unbounded_tiny_qp.h:69-74
)

obj_errors <- function(e, sol, d) {
    x <- sol$x; y <- sol$y; b <- d$data$b; c <- d$data$obj
    switch(e$form,
        info = c(sol$info$pobj - e$opt, sol$info$dobj - e$opt),
        cxby = c(sum(c * x) - e$opt, -sum(b * y) - e$opt),
        quad = {
            xpx <- if (is.null(d$data$P)) 0 else sum(sym_upper_mult(d$data$P, x) * x)
            c(0.5 * xpx + sum(c * x) - e$opt, -0.5 * xpx - sum(b * y) - e$opt)
        },
        planted = c(sum(c * x) - d$upstream$planted[["cx"]], -sum(b * y) + d$upstream$planted[["by"]]))
}

## ------------------------------------------------------------------ replay
rds_files <- sort(list.files("problem_results", pattern = "\\.RDS$"))
res <- list()
for (fname in rds_files) {
    nm <- sub("\\.RDS$", "", fname)
    d <- readRDS(file.path("problem_results", fname))
    e <- expected[[nm]]
    expect_false(is.null(e), info = paste(fname, "has no upstream assertions listed"))
    if (is.null(e)) next
    expect_identical(d$upstream$exitflag, e$flag, info = paste(fname, "saved exit flag"))
    replay <- function() scs(A = d$data$A, P = d$data$P, b = d$data$b, obj = d$data$obj,
                             cone = d$cone, initial = d$initial, control = d$settings)
    if (isTRUE(e$r_error)) {
        ## upstream asserts SCS_FAILED; this package rejects a negative lookback in R first
        expect_error(replay(), pattern = "acceleration_lookback", info = fname)
        next
    }
    sol <- tryCatch(replay(), error = function(err) err)
    expect_false(inherits(sol, "error"), info = fname)
    if (inherits(sol, "error")) next
    res[[nm]] <- sol
    expect_identical(sol$info$status_val, e$flag, info = paste(fname, "exit flag"))
    if (!is.null(e$opt) || e$form == "planted") {
        err <- obj_errors(e, sol, d)
        expect_true(abs(err[1L]) < e$obj_tol && abs(err[2L]) < e$obj_tol,
                    info = paste(fname, "objective error", err[1L], err[2L]))
    }
    if (e$verify) verify_solution_correct(d, sol, e$flag, fname)
    if (!is.null(e$extra)) e$extra(sol, d, res)
}
expect_identical(sort(names(expected)), sort(sub("\\.RDS$", "", rds_files)))
