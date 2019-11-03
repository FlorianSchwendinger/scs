context("PSD - Test")
A <- simple_triplet_matrix(
    i = c(1L, 2L, 3L, 4L, 5L, 7L, 8L, 9L, 1L, 2L, 3L, 5L, 6L, 7L, 8L, 9L, 1L,
          2L, 3L, 4L, 5L, 6L, 7L, 8L, 9L),
    j = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L,
          3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L),
    v = c(-7, -15.556349186104, 3, -21, -15.556349186104, 10, 11.3137084989848,
          5, 7, -25.4558441227157, 8, 14.142135623731, 22.6274169979695,
          -10, -14.142135623731, 3, -2, -11.3137084989848, 1, -5, 2.82842712474619,
          -24.0416305603426, -6, 11.3137084989848, 6), nrow = 9L, ncol = 3L)

rhs <- c(33, -12.7279220613579, 26, 14, 12.7279220613579, 56.5685424949238, 91, 14.142135623731, 15)
obj <- c(1, -1, 1)
cone <- list(f = 0L, l = 0L, s = c(2, 3))

sol <- scs(A = A, b = rhs, obj = obj, cone = cone)
expect_equal(sol$info$status, "Solved")
expect_equal(max(abs(sol$x - c(-0.367666090041563, 1.89832827158511, -0.887550426343585))), 0, tol = 1e-4)
expect_equal(sol$info$pobj, -3.153545, tol = 1e-5)
expect_equal(sol$info$dobj, -3.153545, tol = 1e-5)
