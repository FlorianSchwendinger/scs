## A tiny test

sol  <- scs(A = matrix(c(0.5, 2), 2, 1), b = c(3, 1), obj = 1,
            cone = list(l = 2), control = list(max_iters = 5000))
expect_equal(sol$info$status, "unbounded")
expect_equal(sol$info$pobj, -Inf)
expect_equal(sol$info$dobj, -Inf)
