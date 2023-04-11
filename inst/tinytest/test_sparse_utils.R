## Some tests for utils

s <- Matrix::sparseMatrix(i = c(2,4,3:5), j= c(4,7:5,5), x = 1:5, dims = c(7,7), symmetric = TRUE)

## Test that make_csc_sym_matrix works on plain matrices

ms  <- as.matrix(s)
cm  <- scs:::make_csc_symm_matrix.matrix(ms)
s2 <- Matrix::sparseMatrix(i = cm[["matind"]], p = cm[["matbeg"]], x = cm[["values"]],
                           symmetric = TRUE, index1 = FALSE)

expect_equal(s, s2)


## Test that make_csc_sym_matrix works on simple triplet matrices

ts  <- as(s, "dgTMatrix")
sl  <- slam::simple_triplet_matrix(i = ts@i, j = ts@j, v = ts@x, nrow = 7, ncol = 7)

tm  <- scs:::make_csc_symm_matrix.simple_triplet_matrix(sl)
s3 <- Matrix::sparseMatrix(i = tm[["matind"]], p = tm[["matbeg"]], x = tm[["values"]], dims = c(7, 7),
                           symmetric = TRUE, index1 = FALSE)
expect_equal(s, s3)
