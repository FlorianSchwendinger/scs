## Some tests for utils

s <- Matrix::sparseMatrix(i = c(2,4,3:5), j= c(4,7:5,5), x = 1:5, dims = c(7,7), symmetric = TRUE)

## Test that make_csc_sym_matrix works on plain matrices

ms  <- as.matrix(s)
cm  <- scs:::make_csc_symm_matrix.matrix(ms)
s2 <- Matrix::sparseMatrix(i = cm[["matind"]], p = cm[["matbeg"]], x = cm[["values"]],
                           symmetric = TRUE, index1 = FALSE)

expect_equal(s, s2)


## Test that make_csc_sym_matrix works on simple triplet matrices.
## dgTMatrix slots are 0-based, simple_triplet_matrix indices are 1-based.

ts  <- methods::as(s, "TsparseMatrix")
sl  <- slam::simple_triplet_matrix(i = ts@i + 1L, j = ts@j + 1L, v = ts@x, nrow = 7, ncol = 7)

tm  <- scs:::make_csc_symm_matrix.simple_triplet_matrix(sl)
s3 <- Matrix::sparseMatrix(i = tm[["matind"]], p = tm[["matbeg"]], x = tm[["values"]], dims = c(7, 7),
                           symmetric = TRUE, index1 = FALSE)
expect_equal(s, s3)

## The same through as.simple_triplet_matrix() on the dense matrix, with the
## triplets deliberately out of column order.
sl2 <- slam::as.simple_triplet_matrix(ms)
o <- rev(seq_along(sl2$v)); sl2$i <- sl2$i[o]; sl2$j <- sl2$j[o]; sl2$v <- sl2$v[o]
tm2 <- scs:::make_csc_symm_matrix.simple_triplet_matrix(sl2)
expect_equal(tm, tm2)

## Test the upper-triangle extraction from CsparseMatrix slots (P given as a
## dsCMatrix storing either triangle, or as a full dgCMatrix).
sU <- Matrix::forceSymmetric(s, uplo = "U")
sL <- Matrix::forceSymmetric(s, uplo = "L")
sG <- methods::as(s, "generalMatrix")
cU <- scs:::make_csc_upper_from_csc(sU@i, sU@p, sU@x, 7L, transpose = FALSE)
cL <- scs:::make_csc_upper_from_csc(sL@i, sL@p, sL@x, 7L, transpose = TRUE)
cG <- scs:::make_csc_upper_from_csc(sG@i, sG@p, sG@x, 7L, transpose = FALSE)
expect_equal(cU, cm)
expect_equal(cL, cm)
expect_equal(cG, cm)
