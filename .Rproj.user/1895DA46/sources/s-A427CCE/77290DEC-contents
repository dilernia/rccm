#' Random Covariance Clustering Model
#'
#' This function implements the Random Covariance Clustering Model (RCCM) for joint estimation of
#' sparse precision matrices belonging to multiple clusters or groups. Optimization is conducted using block
#' coordinate descent.
#' @param x List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}
#' @param lambda1 Non-negative scalar. Induces sparsity in subject-level matrices.
#' @param lambda2 Non-negative scalar. Induces similarity between subject-level matrices and cluster-level matrices.
#' @param lambda3 Non-negative scalar. Induces sparsity in cluster-level matrices.
#' @param nclusts Number of clusters or groups.
#' @param delta Threshold for convergence.
#' @param max.iters Maximum number of iterations for block coordinate descent optimization.
#' @param z0s Vector of length \eqn{K} with initial cluster memberships.
#' @return A list of length 3 containing:
#' \enumerate{
#' \item \eqn{p} x \eqn{p} x \eqn{nclusts} array of \eqn{nclusts} number of estimated cluster-level precision matrices (Omega0).
#' \item \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} number of estimated subject-level precision matrices (Omegas).
#' \item \eqn{nclusts} x \eqn{K} matrix of estimated cluster weights for each subject (weights).
#' }
#'
#' @author
#' Andrew DiLernia
#'
#' @examples
#' # Generate data with 2 clusters with 12 and 10 subjects respectively,
#' # 15 variables for each subject, 100 observations for each variable for each subject,
#' # the groups sharing about 50% of network connections, and 10% of differential connections
#' # within each group
#' myData <- rccSim(G = 2, clustSize = c(12, 10), p = 15, n = 100, overlap = 0.50, rho = 0.10)
#'
#' # Analyze simulated data with RCCM
#' result <- rccm(x = myData$simDat, lambda1 = 10, lambda2 = 50, lambda3 = 2, nclusts = 2, delta = 0.001)
#'
#' @export
rccm <- function(x, lambda1, lambda2, lambda3 = 0, nclusts, delta = 0.001, max.iters = 100, z0s = NULL) {

  # Function for making almost symmetric matrix symmetric
  mkSymm <- function(x) {
    return((x + t(x)) / 2)
  }

  # Inputs:
  K <- length(x)
  G <- nclusts
  p <- dim(x[[1]])[2]

  Sl <- sapply(x, cov, simplify = "array")
  nks <- sapply(x, nrow)

  # Initializing subject-level matrices
  Omegas <- sapply(x, FUN = function(datf) {
    mkSymm(glasso::glasso(cov(datf), rho = 0.001)$wi)}, simplify = "array")

  # Initializing weights using hierarchical clustering based on dissimilarity matrix of
  # Frobenius norm of glasso matrix differences
  distMat <- matrix(NA, nrow = K, ncol = K)
  for (r in 1:K) {
    for (s in 1:K) {
      distMat[r, s] <- norm(Omegas[, , r] - Omegas[, , s], type = "F")
    }
  }

  if (is.null(z0s)) {
    cl0 <- cutree(hclust(d = as.dist(distMat), method = "ward.D"), k = G)
  } else {
    cl0 <- z0s
  }

  wgk <- matrix(NA, nrow = G, ncol = K)

  for (i in 1:G) {
    for (j in 1:K) {
      wgk[i, j] <- ifelse(cl0[j] == i, 1, 0)
    }
  }

  # Initializing cluster-level matrices to be all 0's
  Omega0 <- array(0, dim = c(p, p, G))

  Omegas.old <- array(0, c(p, p, K))
  Omega0.old <- array(0, c(p, p, G))
  counter <- 0

  # Initializing vector for deltas and array for weights across iterations
  deltas <- c()
  wArray <- array(NA, dim = c(G, K, max.iters + 1))
  wArray[, , 1] <- wgk

  # Start BCD algorithm
  while (max(abs(Omega0 - Omega0.old)) > delta |
         max(abs(Omegas - Omegas.old)) > delta | counter < 1) {

    # Exit if exceeds max.iters
    if (counter >= max.iters) {
      paste0("Omegas fail to converge for lambda1 =", lambda1, ", lambda2 =", lambda2, ", lambda3 =", lambda3,
          ", delta =", deltas[counter - 1], "\n")

      # Returning results
      res <- list(Omega0, Omegas, wgk)
      names(res) <- c("Omega0", "Omegas", "weights")
      return(res)
    }

    # record current Omega0 & Omegas
    Omega0.old <- Omega0
    Omegas.old <- Omegas

    # 1st step: Updating pi's
    pigs <- 1 / K * rowSums(wgk)

    # 2nd step: updating cluster-level precision matrices

    # Calculating weighted-sum of subject-level matrices
    inv0 <- array(0, c(p, p, G))
    s0 <- sapply(1:G, function(g) {
      wks <- sapply(1:K, function(k) {
        wgk[g, k] * Omegas[, , k]}, simplify = "array")
      return(apply(X = wks, MARGIN = c(1:2), FUN = sum))
    }, simplify = "array")

    for (g in 1:G) {
      S0 <- s0[, , g] / sum(wgk[g, ])
      penMat <- matrix(lambda3 / (lambda2 * sum(wgk[g, ])), nrow = p, ncol = p)
      diag(penMat) <- 0
      if (counter > 1) {
        invisible(capture.output(Omega0[, , g] <- spcov::spcov(Sigma = Omega0[, , g], S = S0, lambda = penMat,
                                                               tol.outer = delta, step.size = 100)$Sigma))
      } else {
        invisible(capture.output(Omega0[, , g] <- spcov::spcov(Sigma = solve(S0), S = S0, lambda = penMat,
                                                               tol.outer = delta, step.size = 100)$Sigma))
      }
      # Calculating inverse of Omega_g for Omega_k and w_gk updates
      inv0[, , g] <- solve(Omega0[, , g])
    }

    # 2b step: updating weights

    # Weight matrix where each row is for a cluster and each column for a subject
    for (i in 1:G) {
      det0 <- det(Omega0[, , i])
      for (j in 1:K) {
        wgk[i, j] <- log(pigs[i]) - lambda2 / 2 * sum(diag(inv0[, , i] %*% Omegas[, , j])) + (-lambda2 / 2) * (log(1 / lambda2^p) + log(det0))
      }
    }
    wgk <- apply(wgk, 2, function(column) {
      column <- exp(column - max(column))
    return(column / sum(column))})

    # 3rd step: updating subject-level precision matrices

    sk <- sapply(1:K, function(k) {
      # Calculating weighted-sum of cluster-level matrices
      ws <- apply(X = sapply(1:G, function(g) {
        wgk[g, k] * inv0[, , g]}, simplify = "array"), MARGIN = c(1:2), FUN = sum)
      return((nks[k] * Sl[, , k] + lambda2 * ws) / (nks[k] + lambda2 - p - 1))
    }, simplify = "array")

    rhoMat <- sapply(X = 1:K, FUN = function(x) {
      matrix(lambda1 / (nks[x] + lambda2 - p - 1), nrow = p, ncol = p)},
                     simplify = "array")

    for (k in 1:K) {
      diag(rhoMat[, , k]) <- 0
      Omegas[, , k] <- mkSymm(glasso::glasso(sk[, , k], rho = rhoMat[, , k], penalize.diagonal = FALSE)$wi)
    }

    # 4th step: updating weights

    # Weight matrix where each row is for a cluster and each column for a subject
    for (i in 1:G) {
      det0 <- det(Omega0[, , i])
      for (j in 1:K) {
        wgk[i, j] <- log(pigs[i]) - lambda2 / 2 * sum(diag(inv0[, , i] %*% Omegas[, , j])) + (-lambda2 / 2) * (log(1 / lambda2^p) + log(det0))
      }
    }
    wgk <- apply(wgk, 2, function(column) {
      column <- exp(column - max(column));
    return(column / sum(column))})

    if (sum(is.nan(wgk)) > 0) {
      stop()
    }

    # Record BCD iterations
    counter <- counter + 1

    deltas <- c(deltas, max(abs(Omega0 - Omega0.old)), max(abs(Omegas - Omegas.old)))
    wArray[, , counter + 1] <- wgk
  }

  # Returning results
  res <- list(Omega0, Omegas, wgk)
  names(res) <- c("Omega0", "Omegas", "weights")
  return(res)
}
