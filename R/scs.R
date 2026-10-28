#  ---------------------------------------------------------
#  scs
#  ===
#' @title SCS - Splitting Conic Solver
#'
#' @description Solves convex cone programs via operator splitting.
#' @param A a matrix of constraint coefficients: a base \code{matrix}, a
#'        \code{dgCMatrix} from package \pkg{Matrix}, or a
#'        \code{simple_triplet_matrix} from package \pkg{slam}. Other
#'        \pkg{Matrix} classes must be converted by the caller.
#'        \bold{NOTE:} The rows of matrix A have to be ordered according to the
#'        order given in subsection \dQuote{Allowed cone parameters}. For more information see \bold{README}.
#' @param b a numeric vector giving the primal constraints
#' @param obj a numeric vector giving the primal objective
#' @param P a symmetric positive semidefinite matrix: a base \code{matrix},
#'        a \code{dsCMatrix} (either triangle stored) or \code{dgCMatrix}
#'        from package \pkg{Matrix}, or a \code{simple_triplet_matrix}. Only
#'        the upper triangle is used. Default \code{NULL}.
#' @param cone a list giving the cone sizes
#' @param initial a named list (warm start solution) of three elements: \code{x} (length = \code{length(obj)}), \code{y} (length = \code{nrow(A)}), and \code{s} (length = \code{nrow(A)}), default NULL indicating no warm start.
#' @param control a list giving the control parameters. For more information see \bold{README}.
#' @return A list with the solution vectors \code{x}, \code{y}, \code{s} and
#'   a list \code{info} describing the run, see section \dQuote{Return
#'   information}.
#' @useDynLib scs, .registration = TRUE
#' @export scs
#' @details
#'
#' \subsection{Important Note}{ \cr
#' The order of the rows in matrix \eqn{A} has to correspond to the order given in
#' the table \dQuote{Cone Arguments}, which means means rows corresponding to
#' \emph{primal zero cones} should be first, rows corresponding to \emph{non-negative cones} second,
#' rows corresponding to \emph{box cones} third, rows corresponding to \emph{second-order cones} fourth,
#' rows corresponding to \emph{positive semidefinite cones} fifth, rows corresponding to
#' \emph{complex positive semidefinite cones} sixth, rows corresponding to \emph{exponential cones}
#' seventh and rows corresponding to \emph{power cones} at last.
#' }
#'
#' \subsection{SCS can solve}{
#' \enumerate{
#'   \item linear programs (LPs)
#'   \item second-order cone programs (SOCPs)
#'   \item semidefinite programs (SDPs), real and complex
#'   \item exponential cone programs (ECPs)
#'   \item power cone programs (PCPs)
#'   \item problems with any combination of cones, which can be defined by the parameters listed in the subsection \dQuote{Allowed cone parameters}
#' } }
#'
#' \subsection{Allowed \emph{cone} parameters are}{
#' \tabular{rllll}{
#'    \tab \bold{Parameter} \tab \bold{Type} \tab \bold{Length} \tab \bold{Description}                       \cr
#'    \tab \code{z}         \tab integer     \tab \eqn{1}       \tab number of primal zero cones (dual free cones),       \cr
#'    \tab                  \tab             \tab               \tab which corresponds to the primal equality constraints \cr
#'    \tab \code{l}         \tab integer     \tab \eqn{1}       \tab number of linear cones (non-negative cones)          \cr
#'    \tab \code{bsize}     \tab integer     \tab \eqn{1}       \tab size of box cone (number of bounds plus one);        \cr
#'    \tab                  \tab             \tab               \tab defaults to \code{length(bl) + 1} when omitted       \cr
#'    \tab \code{bl}        \tab numeric     \tab \eqn{bsize-1} \tab lower limits for box cone (may be \code{-Inf})       \cr
#'    \tab \code{bu}        \tab numeric     \tab \eqn{bsize-1} \tab upper limits for box cone (may be \code{Inf})        \cr
#'    \tab \code{q}         \tab integer     \tab \eqn{\geq1}   \tab vector of second-order cone sizes                    \cr
#'    \tab \code{s}         \tab integer     \tab \eqn{\geq1}   \tab vector of positive semidefinite cone sizes           \cr
#'    \tab \code{cs}        \tab integer     \tab \eqn{\geq1}   \tab vector of complex positive semidefinite cone sizes   \cr
#'    \tab \code{ep}        \tab integer     \tab \eqn{1}       \tab number of primal exponential cones                   \cr
#'    \tab \code{ed}        \tab integer     \tab \eqn{1}       \tab number of dual exponential cones                     \cr
#'    \tab \code{p}         \tab numeric     \tab \eqn{\geq1}   \tab vector of primal/dual power cone parameters
#' }
#' A real positive semidefinite cone of order \eqn{k} occupies \eqn{k(k+1)/2}
#' rows (lower triangle, column-wise, off-diagonal entries scaled by
#' \eqn{\sqrt{2}}); a complex positive semidefinite cone of order \eqn{k}
#' occupies \eqn{k^2} rows (real and imaginary parts of the lower triangle,
#' same scaling). The experimental spectral cones of the SCS C library
#' (log-determinant, nuclear norm, \eqn{\ell_1}, sum of largest eigenvalues)
#' are not compiled into this package. The redundant length fields
#' \code{qsize}, \code{ssize}, \code{cssize} and \code{psize} are accepted
#' for backward compatibility and must equal the length of the vector they
#' describe.
#' }
#'
#' \subsection{Return information}{
#' \code{info} is a list with the elements \code{iter}, \code{status},
#' \code{lin_sys_solver}, \code{status_val}, \code{scale_updates},
#' \code{pobj}, \code{dobj}, \code{res_pri}, \code{res_dual}, \code{gap},
#' \code{res_infeas}, \code{res_unbdd_a}, \code{res_unbdd_p},
#' \code{setup_time}, \code{solve_time} (milliseconds), \code{scale},
#' \code{comp_slack}, \code{rejected_accel_steps},
#' \code{accepted_accel_steps}, \code{aa_stats}, \code{lin_sys_time},
#' \code{cone_time} and \code{accel_time} (milliseconds).
#' \code{status_val} is the SCS exit flag: \code{1} solved, \code{2} solved
#' inaccurately, \code{-1} unbounded, \code{-2} infeasible, \code{-4} failed,
#' \code{-5} interrupted, \code{-6} unbounded inaccurately, \code{-7}
#' infeasible inaccurately. \code{aa_stats} is a list of Anderson
#' acceleration diagnostics: \code{iter}, \code{n_accept},
#' \code{n_reject_lapack}, \code{n_reject_rank0}, \code{n_reject_nonfinite},
#' \code{n_reject_weight_cap}, \code{n_safeguard_reject}, \code{last_rank},
#' \code{last_aa_norm} (\code{NaN} if no acceleration step was attempted)
#' and \code{last_regularization}.
#' }
#'
#' @examples
#' A <- matrix(c(1, 1), ncol=1)
#' b <- c(1, 1)
#' obj <- 1
#' cone <- list(z = 2)
#' control <- list(eps_rel = 1e-3, eps_abs = 1e-3, max_iters = 50)
#' sol <- scs(A = A, b = b, obj = obj, cone = cone, control = control)
#' sol
#  ---------------------------------------------------------
scs <- function(A, b, obj, P = NULL, cone, initial = NULL, control = scs_control()) {

    control <- update_default_controls(control)
    cone <- check_cone(cone)

    ## Objects from package Matrix are handled through their slots only. No
    ## function of Matrix is called (the package does not import Matrix and
    ## must not depend on its coercion API), and dimensions come from @Dim,
    ## never from dim()/NCOL()/NROW(): for a simple_triplet_matrix (a list)
    ## those answer 1, and for an S4 object deserialized before the Matrix
    ## namespace is loaded the first dim() call returns NULL (R bug PR#18207).
    mcls <- matrix_pkg_class(A)
    if (!is.na(mcls)) {
        if (mcls != "dgCMatrix")
            stop("'A' of class '", mcls, "' is not supported; convert it to a ",
                 "'dgCMatrix', e.g. as(as(A, \"generalMatrix\"), \"CsparseMatrix\")")
        n_constraints <- A@Dim[1L]
        n_variables <- A@Dim[2L]
        Ai <- A@i
        Ap <- A@p
        Ax <- A@x
    } else {
        if (inherits(A, "simple_triplet_matrix")) {
            n_constraints <- A$nrow
            n_variables <- A$ncol
        } else {
            A <- as.matrix(A)
            n_constraints <- nrow(A)
            n_variables <- ncol(A)
        }
        csc <- make_csc_matrix(A)
        Ai <- csc[["matind"]]
        Ap <- csc[["matbeg"]]
        Ax <- csc[["values"]]
    }
    if (n_variables < 1L || n_constraints < 1L)
        stop("'A' must have at least one row and one column")
    if (length(b) != n_constraints)
        stop("'b' must have length nrow(A) = ", n_constraints)
    if (length(obj) != n_variables)
        stop("'obj' must have length ncol(A) = ", n_variables)

    data <- list(m = n_constraints, n = n_variables, c = as.double(obj),
                 Ai = Ai, Ap = Ap, Ax = Ax, b = as.double(b))

    if (!is.null(P)) {
      ## SCS wants the upper triangle of P in CSC form. A dsCMatrix stores
      ## one triangle (slot uplo says which); a dgCMatrix may hold the full
      ## symmetric matrix or just its upper triangle. Both are reduced to the
      ## upper triangle in base R (make_csc_upper_from_csc).
      mcls <- matrix_pkg_class(P)
      if (!is.na(mcls)) {
        if (!mcls %in% c("dsCMatrix", "dgCMatrix"))
            stop("'P' of class '", mcls, "' is not supported; convert it to a ",
                 "'dsCMatrix' or 'dgCMatrix', e.g. as(P, \"CsparseMatrix\")")
        pdim <- P@Dim
      } else if (inherits(P, "simple_triplet_matrix")) {
        pdim <- c(P$nrow, P$ncol)
      } else {
        P <- as.matrix(P)
        pdim <- dim(P)
      }
      if (pdim[1L] != n_variables || pdim[2L] != n_variables)
          stop("'P' must be a ", n_variables, " x ", n_variables, " matrix")
      if (!is.na(mcls)) {
        csc <- make_csc_upper_from_csc(P@i, P@p, P@x, n_variables,
                                       transpose = (mcls == "dsCMatrix" && P@uplo == "L"))
      } else {
        csc <- make_csc_symm_matrix(P)
      }
      data$Pi <- csc[["matind"]]
      data$Pp <- csc[["matbeg"]]
      data$Px <- csc[["values"]]
    }

    if (!is.null(initial)) {
        if (!is.list(initial) || !all(c("x", "y", "s") %in% names(initial)))
            stop("'initial' must be a list with elements 'x', 'y' and 's'")
        if (length(initial$x) != n_variables)
            stop("'initial$x' must have length ncol(A) = ", n_variables)
        if (length(initial$y) != n_constraints || length(initial$s) != n_constraints)
            stop("'initial$y' and 'initial$s' must have length nrow(A) = ", n_constraints)
        data$initial <- list(x = as.double(initial$x), y = as.double(initial$y),
                             s = as.double(initial$s))
    }

    ret <- .Call("scsr", data, cone, control, PACKAGE = "scs")
    return(ret)
}

#' @title SCS Control Arguments
#' @description Details to the \emph{control} parameters. The defaults are
#'   those of the bundled SCS C library (see \code{\link{scs_version}}),
#'   except \code{verbose}, which is \code{FALSE} here.
#' @param max_iters an integer giving the maximum number of iterations (default is \code{100000L}).
#' @param eps_rel a double specifying relative feasibility tolerance (default \code{1e-4}).
#' @param eps_abs a double specifying absolute feasibility tolerance (default \code{1e-4}).
#' @param eps_infeas a double specifying infeasibility tolerance (primal and dual) (default \code{1e-7}).
#' @param alpha a double giving the (Douglas-Rachford) over-relaxation parameter, allowed values are in (0, 2) (default \code{1.5}).
#' @param rho_x a double giving the primal scale factor (default \code{1e-6}).
#' @param scale a double giving the initial dual scale factor (default \code{0.1}), updated during the solve if \code{adaptive_scale} is \code{TRUE}.
#' @param verbose a logical giving if the progress should be printed (default is \code{FALSE}).
#' @param normalize a logical giving if heuristic data rescaling should be used (default is \code{TRUE}).
#' @param warm_start a logical indicating if a warm_start is provided (default \code{FALSE}; a call to \code{scs} with a non-null \code{initial} argument overrides it to be effectively \code{TRUE}).
#' @param acceleration_lookback a non-negative integer giving how much memory to use for Anderson acceleration (AA). More memory requires more time to compute but can give more reliable steps; \code{0L} disables acceleration (default \code{10L}). The AA type is chosen by \code{acceleration_type_1}; a negative value, which earlier SCS versions read as a request for type-II AA, is an error.
#' @param acceleration_interval an integer; AA is applied every \code{acceleration_interval} iterations (default \code{5L}).
#' @param adaptive_scale a logical indicating whether to heuristically adapt the dual scale through the solve (default \code{TRUE}).
#' @param write_data_filename a string indicating filename to write problem data to (default \code{NULL} indicating no write).
#' @param log_csv_filename a string indicating filename where SCS will write csv logs of various quantities through the solver (default \code{NULL} indicating no logging, as it makes the solver much slower).
#' @param time_limit_secs a double indicating time limit for solve run in seconds; can be fractional (default \code{0.0} indicating no limit).
#' @param adaptive_diag_scale a logical indicating whether to refine the scaling metric per row from the primal residual profile; requires \code{adaptive_scale} and is silently disabled without it (default \code{TRUE}).
#' @param acceleration_type_1 a logical; \code{TRUE} selects type-I AA, \code{FALSE} type-II, which is more stable but typically slower and is usually run with a smaller \code{acceleration_regularization} such as \code{1e-12} (default \code{TRUE}).
#' @param acceleration_regularization a double giving the Tikhonov regularization of the AA least-squares solve. A positive value is scaled by the Frobenius norms of the AA matrices, a negative value pins the absolute regularization at \code{-acceleration_regularization}, and \code{0} disables regularization (default \code{1e-8}).
#' @param acceleration_relaxation a double in \eqn{[0, 2]} giving the AA relaxation factor; \code{1.0} is vanilla AA (default \code{1.0}).
#' @return a list containing the control parameters.
#' @export scs_control
scs_control <- function(max_iters = 100000L, eps_rel = 1e-4, eps_abs = 1e-4, eps_infeas = 1e-7,
                        alpha = 1.5, rho_x = 1e-6, scale = 0.1, verbose = FALSE, normalize = TRUE,
                        warm_start = FALSE, acceleration_lookback = 10L, acceleration_interval = 5L,
                        adaptive_scale = TRUE, write_data_filename = NULL, log_csv_filename = NULL,
                        time_limit_secs = 0.0, adaptive_diag_scale = TRUE,
                        acceleration_type_1 = TRUE, acceleration_regularization = 1e-8,
                        acceleration_relaxation = 1.0) {
    as.list(environment())
}

update_default_controls <- function(control) {
    cntrl <- scs_control()
    unknown <- setdiff(names(control), names(cntrl))
    if (length(unknown))
        warning("ignoring unknown control parameter(s): ",
                paste(sQuote(unknown), collapse = ", "), call. = FALSE)
    nam <- names(control)[names(control) %in% names(cntrl)]
    cntrl[nam] <- control[nam]
    ## SCS 3.3 rejects a negative lookback (it would fail the solve with status
    ## -4 and NaN results); say what replaced the old convention.
    lb <- cntrl[["acceleration_lookback"]]
    if (length(lb) == 1L && !is.na(lb) && lb < 0)
        stop("'acceleration_lookback' must be non-negative; to select type-II Anderson ",
             "acceleration set 'acceleration_type_1 = FALSE' instead of a negative lookback",
             call. = FALSE)
    cntrl
}

## Class name of an S4 object from package Matrix, NA otherwise. A pure
## attribute test: no S4 dispatch, so the Matrix namespace is never loaded by
## this package (see the note in scs()).
matrix_pkg_class <- function(x) {
    if (isS4(x) && identical(attr(class(x), "package"), "Matrix")) class(x)[1L]
    else NA_character_
}

.cone_names <- c("z", "l", "bsize", "bl", "bu", "q", "s", "cs", "ep", "ed", "p")
## Redundant array lengths accepted for backward compatibility (the C struct
## has them; earlier versions of this package ignored them). They are checked
## against the array they describe and then dropped.
.cone_size_names <- c(q = "qsize", s = "ssize", cs = "cssize", p = "psize")
.spectral_cone_names <- c("d", "nuc_m", "nuc_n", "ell1", "sl_n", "sl_k")

check_cone <- function(cone) {
    if (!is.list(cone)) stop("'cone' must be a named list")
    nms <- names(cone)
    if (length(cone) && (is.null(nms) || any(!nzchar(nms))))
        stop("'cone' must be a named list")
    ## Empty entries for cones this package does not handle are dropped, as
    ## 3.2.7 did: CVXPY-style callers always send e.g. 'pnd' (an empty list).
    ## A non-empty one is an error rather than a silently ignored cone.
    unused <- setdiff(nms, c(.cone_names, .cone_size_names))
    cone[unused[lengths(cone[unused]) == 0L]] <- NULL
    nms <- names(cone)
    spectral <- intersect(nms, .spectral_cone_names)
    if (length(spectral))
        stop("spectral cone(s) ", paste(sQuote(spectral), collapse = ", "),
             " are not compiled into this package")
    unknown <- setdiff(nms, c(.cone_names, .cone_size_names))
    if (length(unknown))
        stop("unknown cone parameter(s): ", paste(sQuote(unknown), collapse = ", "))
    for (arr in names(.cone_size_names)) {
        size_name <- .cone_size_names[[arr]]
        if (!is.null(cone[[size_name]])) {
            if (length(cone[[size_name]]) != 1L || cone[[size_name]] != length(cone[[arr]]))
                stop("'", size_name, "' must equal length(", arr, ") = ", length(cone[[arr]]))
            cone[[size_name]] <- NULL
        }
    }
    ## Box cone: 'bsize' counts the bounds plus the scaling variable t. When
    ## it is omitted it is derived from the bounds, so that 'bl'/'bu' alone
    ## define the cone instead of being silently ignored.
    nbl <- length(cone[["bl"]]); nbu <- length(cone[["bu"]])
    if (nbl != nbu)
        stop("box cone: 'bl' and 'bu' must have the same length (", nbl, " vs ", nbu, ")")
    bsize <- cone[["bsize"]]
    if (is.null(bsize)) {
        if (nbl > 0) cone[["bsize"]] <- nbl + 1L
    } else {
        if (length(bsize) != 1L || is.na(bsize) || bsize < 0)
            stop("box cone: 'bsize' must be a single non-negative integer")
        if (bsize > 0 && nbl != bsize - 1)
            stop("box cone: 'bl' and 'bu' must have length bsize - 1 = ", bsize - 1)
        if (bsize == 0 && nbl > 0)
            stop("box cone: 'bl'/'bu' given but 'bsize' is 0; omit 'bsize' or set it to length(bl) + 1")
    }
    cone
}

#' @title Version of the Bundled SCS Library
#' @description Returns the version string of the SCS C library compiled
#'   into this package.
#' @return a character string, for example \code{"3.3.1"}.
#' @examples
#' scs_version()
#' @export scs_version
scs_version <- function() {
    .Call("scsr_version", PACKAGE = "scs")
}
