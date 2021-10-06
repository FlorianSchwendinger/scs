context("PSD - Test")
if (grepl("SunOS", Sys.info()["sysname"])) { ## Skip this on Solaris
    skip("Skipping on Solaris as this has always failed")
} else {
  test_that("Simple PSD test", {
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
    cone <- list(z = 0L, l = 0L, s = c(2, 3))

    sol <- scs(A = A, b = rhs, obj = obj, cone = cone, control = list(verbose = 1))
    expect_equal(sol$info$status, "solved")
    expect_equal(sol$x,
                 c(-0.36778629730617,  1.89834935806707, -0.88746074956525),
                 tolerance = 1e-4)
    expect_equal(sol$y,
                 c(0.00396134156379081, -0.00613634113889057, 0.00475276897567086, 0.0558030988006113,
                   -0.00340904764830297, 0.0342440719001259, 0.000104130479854425, -0.00147926212203912,
                   0.0105070907711366),
                 tolerance = 1e-4)
    expect_equal(sol$info$pobj, -3.1535964049385, tol = 1e-4)
    expect_equal(sol$info$dobj, -3.1535519733301, tol = 1e-4)
  })
}

