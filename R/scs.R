#  ---------------------------------------------------------
#  scs
#  ===
#' @title SCS - Splitting Conic Solver
#'
#' @description Solves convex cone programs via operator splitting.
#' @param A a matrix of constraint coefficients.
#'        \bold{NOTE:} The rows of matrix A have to be ordered according to the
#'        order given in subsection \dQuote{Allowed cone parameters}. For more information see \bold{README}.
#' @param b a numeric vector giving the primal constraints
#' @param obj a numeric vector giving the primal objective
#' @param P a symmetric positive semidefinite matrix, default \code{NULL}
#' @param cone a list giving the cone sizes
#' @param initial a named list (warm start solution) of three elements: \code{x} (length = \code{length(obj)}), \code{y} (length = \code{nrow(A)}), and \code{s} (length = \code{nrow(A)}), default NULL indicating no warm start.
#' @param control a list giving the control parameters. For more information see \bold{README}.
#' @return list of solution vectors x, y, s and information about run
#' @useDynLib scs, .registration = TRUE
#' @export scs
#' @details
#'
#' \subsection{Important Note}{ \cr
#' The order of the rows in matrix \eqn{A} has to correspond to the order given in
#' the table \dQuote{Cone Arguments}, which means means rows corresponding to
#' \emph{primal zero cones} should be first, rows corresponding to \emph{non-negative cones} second,
#' rows corresponding to \emph{second-order cone} third, rows corresponding to \emph{positive semidefinite cones} fourth,
#' rows corresponding to \emph{exponential cones} fifth and rows corresponding to \emph{power cones} at last.
#' }
#'
#' \subsection{SCS can solve}{
#' \enumerate{
#'   \item linear programs (LPs)
#'   \item second-order cone programs (SOCPs)
#'   \item semidefinite programs (SDPs)
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
#'    \tab \code{bsize}     \tab integer     \tab \eqn{1}       \tab size of box cone                                     \cr
#'    \tab \code{bl}        \tab numeric     \tab \eqn{bsize-1} \tab lower limit for box cone                             \cr
#'    \tab \code{bu}        \tab numeric     \tab \eqn{bsize-1} \tab upper limit for box cone                             \cr
#'    \tab \code{q}         \tab integer     \tab \eqn{\geq1}   \tab vector of second-order cone sizes                    \cr
#'    \tab \code{s}         \tab integer     \tab \eqn{\geq1}   \tab vector of positive semidefinite cone sizes           \cr
#'    \tab \code{ep}        \tab integer     \tab \eqn{1}       \tab number of primal exponential cones                   \cr
#'    \tab \code{ed}        \tab integer     \tab \eqn{1}       \tab number of dual exponential cones                     \cr
#'    \tab \code{p}         \tab numeric     \tab \eqn{\geq1}   \tab vector of primal/dual power cone parameters
#' } }
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

    ## TBD check box cone parameters, bsize > 0  & bl, bu have lengths bsize - 1

    n_variables <- NCOL(A)
    n_constraints <- NROW(A)

    if ( inherits(A, "dgCMatrix") ) {
        Ai <- A@i
        Ap <- A@p
        Ax <- A@x
    } else {
        csc <- make_csc_matrix(A)
        Ai <- csc[["matind"]]
        Ap <- csc[["matbeg"]]
        Ax <- csc[["values"]]
    }

    data <- list(m = n_constraints, n = n_variables, c = obj,
                 Ai = Ai, Ap = Ap, Ax = Ax, b = b)

    if (!is.null(P)) {
      if (inherits(P, "dsCMatrix") ) {
        data$Pi <- P@i
        data$Pp <- P@p
        data$Px <- P@x
      } else {
        csc  <- make_csc_symm_matrix(P)
        data$Pi <- csc[["matind"]]
        data$Pp <- csc[["matbeg"]]
        data$Px <- csc[["values"]]
      }
    }

    data$initial <- initial

    ret <- .Call("scsr", data, cone, control, PACKAGE = "scs")
    return(ret)
}

#' @title SCS Control Arguments
#' @description Details to the \emph{control} parameters.
#' @param max_iters an integer giving the maximum number of iterations (default is \code{100000L}).
#' @param eps_rel a double specifying relative feasibility tolerance (default \code{1e-4}).
#' @param eps_abs a double specifying absolute feasibility tolerance (default \code{1e-4}).
#' @param eps_infeas a double specifying infeasibility tolerance (primal and dual) (default \code{1e-7}).
#' @param alpha a double giving the (Douglas-Rachford) over-relaxation parameter, allowed values are in (0, 2) (default \code{1.5}).
#' @param rho_x a double giving the momentum of x term (default os \code{1e-6}).
#' @param scale a double giving the factor (default is \code{1.0}) by which the data is rescaled (only used if normalize is \code{TRUE}).
#' @param verbose a logical giving if the progress should be printed (default is \code{FALSE}).
#' @param normalize a logical giving if heuristic data rescaling should be used (default is \code{TRUE}).
#' @param warm_start a logical indicating if a warm_start is provided (default \code{FALSE}, but a call to \code{scs} with a non-null \code{initial} argument overrides it to be effectively \code{TRUE})
#' @param acceleration_lookback an integer indicating How much memory to use for Anderson acceleration. More memory requires more time to compute but can give more reliable steps (default \code{0L}, disabling it).
#' @param acceleration_interval an integer specifying the number of iterations for which Anderson acceleration is run (default \code{1L}).
#' @param adaptive_scale a logical indicating whether to heuristically adapt dual through the solve (default \code{TRUE}).
#' @param write_data_filename a string indicating filename to write problem data to (default \code{NULL} indicating no write).
#' @param log_csv_filename a string indicating filename where SCS will write csv logs of various quantities through the solver (default \code{NULL} indicating no logging, as it makes the solver much slower).
#' @param time_limit_secs a double indicating time limit for solve run in seconds; can be fractional (default \code{0.0} indicating no limit).
#' @return a list containing the control parameters.
#' @export scs_control
scs_control <- function(max_iters = 100000L, eps_rel = 1e-4, eps_abs = 1e-4, eps_infeas = 1e-7,
                        alpha = 1.5, rho_x = 1e-6, scale = 0.1, verbose = FALSE, normalize = TRUE,
                        warm_start = FALSE, acceleration_lookback = 0L, acceleration_interval = 1L,
                        adaptive_scale = TRUE, write_data_filename = NULL, log_csv_filename = NULL,
                        time_limit_secs = 0.0) {
    as.list(environment())
}

update_default_controls <- function(control) {
    cntrl <- scs_control()
    nam <- names(control)[names(control) %in% names(cntrl)]
    cntrl[nam] <- control[nam]
    cntrl
}
