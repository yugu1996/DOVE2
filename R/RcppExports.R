# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

BS <- function(t, knots) {
    .Call(`_DOVE2_BS`, t, knots)
}

newCox <- function(beta, gamma, time, delta, X, W, S, knots, threshold = 10^(-4), maxit = 500L) {
    .Call(`_DOVE2_newCox`, beta, gamma, time, delta, X, W, S, knots, threshold, maxit)
}

newCox_noX <- function(gamma, time, delta, W, S, knots, threshold = 10^(-4), maxit = 500L) {
    .Call(`_DOVE2_newCox_noX`, gamma, time, delta, W, S, knots, threshold, maxit)
}

