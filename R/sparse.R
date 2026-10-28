##
## copied from the Rsymphony package from
## Kurt Hornik and Stefan Theussl and Reinhard Harter
##
## Simple functions for converting "matrix" type objects into the
## sparse "column major order" (CSC, modulo offsets) format used by
## SYMPHONY.

## matind: vector of the row indices corresponding to each entry of
##   value
## values: vector of the values of nonzero entries of the constraint
##   matrix in column order.

make_csc_matrix <- function(x) UseMethod("make_csc_matrix")

make_csc_matrix.matrix <- function(x) {
    if(!is.matrix(x))
        stop("Argument 'x' must be a matrix.")
    if (anyNA(x))
        stop("the matrix must not contain NA values")

    ## which(arr.ind = TRUE) returns the entries in column-major order,
    ## which is exactly CSC order.
    ind <- which(x != 0, arr.ind = TRUE)
    list(matbeg = c(0L, cumsum(tabulate(ind[, 2L], ncol(x)))),
         matind = ind[, 1] - 1L,
         values = x[ind])
}

make_csc_matrix.simple_triplet_matrix <- function(x) {
    if(!inherits(x, "simple_triplet_matrix"))
        stop("Argument 'x' must be of class 'simple_triplet_matrix'.")

    ## The matrix method assumes that indices for non-zero entries are
    ## in row-major order, but the simple_triplet_matrix() constructor
    ## currently does not canonicalize accordingly ...
    ind <- order(x$j, x$i)
    list(matbeg = c(0L, cumsum(tabulate(x$j[ind], x$ncol))),
         matind = x$i[ind] - 1L,
         values = x$v[ind])
}

## Upper triangle of a symmetric matrix in CSC form, for SCS's P.
## Added by @bnaras; rewritten for the 3.3.1 release: the triplet method
## added 1 to indices that are already 1-based and neither method reordered
## the values along with the indices.

## 1-based triplets (i, j, v) -> CSC with 0-based row indices, entries sorted
## by column then row, values permuted alongside.
csc_from_triplets <- function(i, j, v, ncol) {
    ord <- order(j, i)
    list(matbeg = c(0L, cumsum(tabulate(j[ord], ncol))),
         matind = as.integer(i[ord]) - 1L,
         values = v[ord])
}

make_csc_symm_matrix <- function(m) UseMethod("make_csc_symm_matrix")

make_csc_symm_matrix.matrix  <- function(m) {
    if (anyNA(m))
        stop("'P' must not contain NA values")
    ind <- which(m != 0, arr.ind = TRUE)
    ind <- ind[ind[, 1L] <= ind[, 2L], , drop = FALSE] ## keep upper part only
    csc_from_triplets(ind[, 1L], ind[, 2L], m[ind], ncol(m))
}

make_csc_symm_matrix.simple_triplet_matrix  <- function(m) {
    keep <- m$i <= m$j  ## slam indices are 1-based already
    csc_from_triplets(m$i[keep], m$j[keep], m$v[keep], m$ncol)
}

## From the slots of a CsparseMatrix (0-based i, p; values x; n columns) to
## the upper-triangle CSC that SCS expects. With transpose = TRUE the stored
## entries are the lower triangle of a symmetric matrix (dsCMatrix with
## uplo = "L") and are mirrored; otherwise entries below the diagonal are
## dropped (a full symmetric dgCMatrix, or one already upper triangular).
make_csc_upper_from_csc <- function(i, p, x, n, transpose = FALSE) {
    j <- rep.int(seq_len(n), diff(p))  ## 1-based column of every entry
    i <- i + 1L
    if (transpose) {
        tmp <- i; i <- j; j <- tmp
    } else {
        keep <- i <= j
        i <- i[keep]; j <- j[keep]; x <- x[keep]
    }
    csc_from_triplets(i, j, x, n)
}

