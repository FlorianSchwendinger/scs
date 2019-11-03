context("Tiny test")

sol  <- scs(matrix(c(0.5,2),2,1), c(3, 1), 1, list(l=2), list(max_iters=5000))
expect_equal(sol$info$status, "Unbounded")
expect_equal(sol$info$pobj, -Inf)
expect_equal(sol$info$dobj, -Inf)
