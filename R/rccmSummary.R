#' Probability Density Function for Wishart Distribution
#'
#' This function provides the probability density function for
#' the Wishart distribution.
#' @param x \eqn{p} x \eqn{p} positive definite matrix.
#' @param M \eqn{p} x \eqn{p} mean matrix. Note that \eqn{M = nu*V} where \eqn{V} is the scale matrix.
#' @param nu Degrees of freedom.
#' @param logged Logical. If TRUE, probability given on log scale.
#' @return The probability density function value.
#'
#' @export
dwishart <- function(x, M, nu, logged = FALSE) {
  x <- (x + t(x)) / 2
  M <- (M + t(M)) / 2
  p <- nrow(x)
  lnumr <- (nu - p - 1) / 2 * log(det(x)) - nu / 2 * sum(diag(solve(M) * x))
  ldenom <- (nu * p / 2) * log(2) + (nu / 2) * log(det(1 / nu * M)) + (p * (p - 1) / 4) * log(pi) + sum(sapply(1:p,
                                                                                         FUN = function(j) {
                                                                                           lgamma(nu / 2 + (1 - j) / 2)}))
  if (logged) {
    return(lnumr - ldenom)
  } else {
    return(exp(lnumr - ldenom))
  }
}

#' Adjacency Matrix
#'
#' This function calculates an adjacency matrix for the
#' matrix \eqn{mat} based on the absolute value threshold of \eqn{thresh}.
#' @param mat Numeric matrix.
#' @param thresh Threshold for absolute value of entries.
#' @return An adjacency matrix containing 1's and 0's.
#'
#' @export
adj <- function(mat, thresh = 0.001) {
  return((abs(mat) > thresh) + 0)
}

#' Z to Adjacency Matrix
#'
#' This function calculates an adjacency matrix based on
#' an integer vector of cluster memberships.
#' @param z Integer vector of cluster memberships.
#' @return An adjacency matrix containing 1's and 0's.
#'
#' @examples
#' # Calculate adjacency matrix for clustering
#' zToA(c(rep(1, 2), rep(2, 2), rep(3, 2)))
#'
#' @export
zToA <- function(z) {
  K <- length(z)
  A <- matrix(0, nrow = K, ncol = K)
  for (r in 1:K) {
    for (s in 1:K) {
      A[r, s] <- ifelse(z[r] != 0 & z[s] != 0, as.integer(z[r] == z[s]), 0)
    }
  }
  return(A)
}

#' Rand Index
#'
#' This function calculates the rand index describing
#' the amount of agreement between two integer vectors
#'  of cluster memberships.
#' @param x First integer vector of cluster memberships.
#' @param y Second integer vector of cluster memberships.
#' @return The rand index value, bounded between 0 and 1.
#'
randCalc <- function(x, y) {
  Ahat <- zToA(x)[lower.tri(zToA(x), diag = FALSE)]
A0 <- zToA(y)[lower.tri(zToA(y), diag = FALSE)]
return((sum((Ahat - A0) == 2) + sum((Ahat - A0) == 0)) / choose(n = length(x), 2))
}

#' Modified AIC and BIC
#'
#' This function calculates modified AIC and BIC values
#' for the random covariance clustering model (RCCM)
#' @param x List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}.
#' @param omegaks \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} number of estimated subject-level precision matrices.
#' @param omega0s \eqn{p} x \eqn{p} x \eqn{nclusts} array of \eqn{nclusts} number of estimated cluster-level precision matrices.
#' @param lambda2 Scalar tuning parameter value used for RCCM.
#' @param ws \eqn{nclusts} x \eqn{K} matrix of estimated cluster weights for each subject.
#' @return Numeric vector of length 2 containing the modified AIC and BIC values.
#'
#' @examples
#' # Generate data
#' set.seed(1994)
#' myData <- rccSim(G = 2, clustSize = 10, p = 10, n = 100, overlap = 0.50, rho = 0.10)
#'
#' # Analyze with RCCM
#' resultRccm <- rccm(x = myData$simDat, lambda1 = 20,
#' lambda2 = 325, lambda3 = 0.01, nclusts = 2)
#'
#' # Calculate modified AIC and BIC
#' aicbicm(x = myData$simDat, omegaks = resultRccm$Omegas,
#' omega0s = resultRccm$Omega0, lambda2 = 325, ws = resultRccm$weights)
#'
#' @export
aicbicm <- function(x, omegaks, omega0s, lambda2, ws) {

  K <- dim(omegaks)[3]
  nks <- sapply(x, FUN = nrow)
  G <- dim(omega0s)[3]
  pigs <- rowSums(ws)

  dfks <- sapply(X = 1:K, FUN = function(x) {
    sum(adj(omegaks[, , x])[lower.tri(omegaks[, , x])])})
  dfgs <- sapply(X = 1:G, FUN = function(x) {
    sum(adj(omega0s[, , x])[lower.tri(omega0s[, , x])])})

  Sk <- lapply(X = 1:K, FUN = function(k){cov(x[[k]])*nks[k]})

  nll <- sapply(X = 1:K, FUN = function(x) {
    0.50*sum(diag(Sk[[x]] %*% omegaks[, , x])) - 0.50*log(det(omegaks[, , x])) - 2 * sum(sapply(1:G, FUN = function(g) {
      ws[g, x] * (log(pigs[g] + dwishart(x = omegaks[, , x], M = omega0s[, , g], nu = lambda2, logged = TRUE)))
      }))
  }) * nks
  bicm <- sum(nll) + sum(log(nks) * dfks) + sum(log(K * pigs) * dfgs)
  aicm <- sum(nll) + sum(2 * dfks) + sum(2 * dfgs)
  return(setNames(c(aicm, bicm), c("AICm", "BICm")))
}

#' AIC and BIC
#'
#' This function calculates AIC and BIC values
#' for the random covariance clustering model (RCCM)
#' @param x List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}.
#' @param omegaks \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} number of estimated subject-level precision matrices.
#' @return Numeric vector of length 2 containing the AIC and BIC values.
#'
#' @examples
#' # Generate data
#' set.seed(1994)
#' myData <- rccSim(G = 2, clustSize = 10, p = 10, n = 100, overlap = 0.50, rho = 0.10)
#'
#' # Analyze with RCCM
#' resultRccm <- rccm(x = myData$simDat, lambda1 = 20,
#' lambda2 = 325, lambda3 = 0.01, nclusts = 2)
#'
#' # Calculate AIC and BIC
#' aicbic(x = myData$simDat, omegaks = resultRccm$Omegas)
#'
#' @export
aicbic <- function(x, omegaks) {

  K <- dim(omegaks)[3]
  nks <- sapply(x, FUN = nrow)

  dfks <- sapply(X = 1:K, FUN = function(x) {
    sum(adj(omegaks[, , x])[lower.tri(omegaks[, , x])])
    })

  Sk <- lapply(X = 1:K, FUN = function(k){cov(x[[k]])*nks[k]})

  nll <- sapply(X = 1:K, FUN = function(x) {
    0.50*sum(diag(Sk[[x]] %*% omegaks[, , x])) - 0.50*log(det(omegaks[, , x]))
  }) * nks
  bic <- sum(nll) + sum(log(nks) * dfks)
  aic <- sum(nll) + sum(2 * dfks)
  return(setNames(c(aic, bic), c("AIC", "BIC")))
}
