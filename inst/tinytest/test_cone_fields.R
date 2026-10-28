## Cone fields this package does not handle: empty ones are dropped (CVXPY-style
## callers always send 'pnd'), non-empty ones are an error.

A <- matrix(c(1, 1), ncol = 1)
b <- c(1, 1)
obj <- -1
ctrl <- list(eps_abs = 1e-9, eps_rel = 1e-9)

ref <- scs(A = A, b = b, obj = obj, cone = list(l = 2), control = ctrl)
expect_equal(ref$info$status, "solved")
expect_equal(ref$x, 1, tolerance = 1e-8)

sol <- scs(A = A, b = b, obj = obj, cone = list(l = 2, pnd = list()), control = ctrl)
expect_identical(sol$x, ref$x)

sol <- scs(A = A, b = b, obj = obj, cone = list(l = 2, d = integer(0)), control = ctrl)
expect_identical(sol$x, ref$x)

expect_error(scs(A = A, b = b, obj = obj, cone = list(l = 2, pnd = list(c(0.5, 0.5))), control = ctrl),
             "unknown cone parameter")
expect_error(scs(A = A, b = b, obj = obj, cone = list(l = 2, d = 2L), control = ctrl),
             "not compiled into this package")
