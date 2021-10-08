context("Simple Example")

test_that("A simple example", {
  A <- matrix(c(1), 2, 1)
  b  <-  c(1,1)
  c  <-  c(1)
  cone  <-  list(z = 2)
  params  <-  list(eps_abs = 1e-3, eps_rel = 1e-3, max_iters = 50)
  sol <- scs(A = A, b = b, obj = c, cone = cone, control = params)

  expect_equal(sol$info$status, "solved")
  expect_equal(sol$x, 1)
  expect_equal(sol$y, c(-0.5, -0.5))
  expect_equal(sol$info$pobj, 1)
})
