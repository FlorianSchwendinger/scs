## An all-zero P (no stored entries) is a valid quadratic term, the same
## problem as P = NULL. CVXR sends one for LPs posed with a quadratic objective.

A <- matrix(c(1, 1), ncol = 1)
b <- c(1, 1)
obj <- -1
cone <- list(l = 2)
ctrl <- list(eps_abs = 1e-9, eps_rel = 1e-9)

ref <- scs(A = A, b = b, obj = obj, cone = cone, control = ctrl)
expect_equal(ref$info$status, "solved")

sol <- scs(A = A, b = b, obj = obj, P = matrix(0, 1, 1), cone = cone, control = ctrl)
expect_equal(sol$info$status, "solved")
expect_identical(sol$x, ref$x)

if (requireNamespace("Matrix", quietly = TRUE)) {
    P0 <- new("dsCMatrix", i = integer(0), p = c(0L, 0L), x = numeric(0),
              Dim = c(1L, 1L), uplo = "U")
    sol <- scs(A = A, b = b, obj = obj, P = P0, cone = cone, control = ctrl)
    expect_equal(sol$info$status, "solved")
    expect_identical(sol$x, ref$x)
}
