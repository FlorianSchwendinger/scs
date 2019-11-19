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
#' @param cone a list giving the cone sizes
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
#'    \tab \code{f}         \tab integer     \tab \eqn{1}       \tab number of primal zero cones (dual free cones),       \cr
#'    \tab                  \tab             \tab               \tab which corresponds to the primal equality constraints \cr
#'    \tab \code{l}         \tab integer     \tab \eqn{1}       \tab number of linear cones (non-negative cones)          \cr
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
#' cone <- list(f = 2)
#' control <- list(eps = 1e-3, max_iters = 50)
#' sol <- scs(A, b, obj, cone, control)
#' sol
#  max_iters=2500L, normalize=TRUE, verbose=TRUE,
#  cg_rate=2.0, scale=5.0, rho_x=1e-03, alpha=1.5, eps=1e-3
#  ---------------------------------------------------------
scs <- function(A, b, obj, cone, control = scs_control()) {

    control <- update_default_controls(control)

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

    ret <- .Call("scsr", data, cone, control, PACKAGE = "scs")
    return(ret)
}

#' @title SCS Control Arguments
#' @description Details to the \emph{control} parameters.
#' @param max_iters an integer giving the maximum number of iterations (default is \code{5000L}).
#' @param normalize a logical giving if heuristic data rescaling should be used (default is \code{TRUE}).
#' @param verbose a logical giving if the progress should be printed (default is \code{FALSE}).
#' @param cg_rate a double giving the rate at which the CG tolerance for the indirect method is tightened (higher is tighter, default is \code{2.0}).
#' @param scale a double giving the factor (default is \code{1.0}) by which the data is rescaled (only used if normalize is \code{TRUE}).
#' @param rho_x a double giving the momentum of x term (default os \code{1e-3}).
#' @param alpha a double giving the over-relaxation parameter, allowed values are in (0, 2) (default if \code{1.5}).
#' @param eps a double giving the convergence tolerance (default is \code{1e-5}).
#' @param acceleration_lookback an integer indicating the number of iterations to look back for Anderson acceleration (default is \code{10L}).
#' @return a list containing the control parameters.
#' @export scs_control
scs_control <- function(max_iters = 5000L, normalize = TRUE, verbose = FALSE,
                        cg_rate = 2.0, scale = 1.0, rho_x = 1e-03, alpha = 1.5, eps = 1e-5,
                        acceleration_lookback = 10L) {
    as.list(environment())
}

update_default_controls <- function(control) {
    cntrl <- scs_control()
    nam <- names(control)[names(control) %in% names(cntrl)]
    cntrl[nam] <- control[nam]
    cntrl
}
