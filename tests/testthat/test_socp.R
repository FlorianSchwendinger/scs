context("SOCP Test")

test_that("An SOCP Test", {
  A <- simple_triplet_matrix(i = 1:4, j = c(1L, 1L, 2L, 3L), v = c(1, -1, -1, -1), nrow = 4L, ncol = 3L)
  rhs <- c(1.4142135623731, 0, 0, 0)
  obj <- c(1, 1, 1)
  cone <- list(z = 1L, l = 0L, q = 3L)

  sol <- scs(A = A, b = rhs, obj = obj, cone = cone)
  expect_equal(sol$info$status, "solved")
  expect_equal(sol$x, c(sqrt(2), -1, -1), tolerance = 1e-7)
  expect_equal(sol$info$pobj, -0.5857864, tol = 1e-5)
  expect_equal(sol$info$dobj, -0.5857864, tol = 1e-5)
})
