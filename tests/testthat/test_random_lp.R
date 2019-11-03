context("Random LP")
set.seed(1)
m  <-  100
n  <- 100
A <- matrix(rnorm(m*n),m,n)
b  <-  rnorm(m)
c  <-  rnorm(n)
cone  <-  list(f = m)
params  <-  list(eps = 1e-4, max_iters = 5000)
sol <- scs(A, b, c, cone, params)
expect_equal(sol$info$status, "Solved")
expect_equal(sol$info$pobj, 13.5089, tol=1e-4)
expect_equal(sol$info$dobj, 13.5089, tol=1e-4)
