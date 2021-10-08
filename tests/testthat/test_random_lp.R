context("Random LP")
test_that("A Random LP Problem", {
  set.seed(1)
  m  <-  100
  n  <- 100
  A <- matrix(rnorm(m*n),m,n)
  b  <-  rnorm(m)
  c  <-  rnorm(n)
  cone  <-  list(z = m)
  params  <-  list(eps_rel = 1e-4, eps_abs = 1e-4, max_iters = 5000)
  sol <- scs(A = A, b = b, obj = c, cone = cone, control = params)
  expect_equal(sol$info$status, "solved")
  expect_equal(sol$info$pobj, 13.5089, tol=1e-4)
  expect_equal(sol$info$dobj, 13.5089, tol=1e-4)
})
