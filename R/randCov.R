#' Random Covariance Model
#'
#' This function implements the Random Covariance Model (RCM) for joint estimation of
#' multiple sparse precision matrices. Optimization is conducted using block
#' coordinate descent.
#' @param x List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}.
#' @param lambda1 Non-negative scalar. Induces sparsity in subject-level matrices.
#' @param lambda2 Non-negative scalar. Induces similarity between subject-level matrices and group-level matrix.
#' @param lambda3 Non-negative scalar. Induces sparsity in group-level matrix.
#' @param delta Threshold for convergence.
#' @param max.iters Maximum number of iterations for block coordinate descent optimization.
#' @return A list of length 2 containing:
#' \enumerate{
#' \item Group-level precision matrix estimate (Omega0).
#' \item \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} subject-level precision matrix estimates (Omegas).
#' }
#' @author
#' Lin Zhang and Andrew DiLernia
#'
#' @examples
#' # Generate data with 5 subjects, 15 variables for each subject,
#' # 100 observations for each variable for each subject,
#' # and 10% of differential connections
#' # within each group
#' myData <- rccSim(G = 1, clustSize = 5, p = 15, n = 100, rho = 0.10)
#'
#' # Analyze simulated data with RCM
#' result <- randCov(x = myData$simDat, lambda1 = 0.30, lambda2 = 0.10, lambda3 = 0.001, delta = 0.001)
#'
#' @export
#'
#' @references
#' Zhang, Lin, Andrew DiLernia, Karina Quevedo, Jazmin Camchong, Kelvin Lim, and Wei Pan.
#' "A Random Covariance Model for Bi-level Graphical Modeling with Application to Resting-state FMRI Data." 2019. https://arxiv.org/pdf/1910.00103.pdf

randCov <- function(x, lambda1, lambda2, lambda3 = 0,
                    delta = 0.001, max.iters = 100) {

  # Sparse covariance optimization function
  spcov_bcd <- function(samp_cov, rho, initial = NULL, lambda2, lambda3) {
    p <- dim(samp_cov)[1]

    if (is.null(initial)) {
      Sigma <- samp_cov + 0.01 * diag(p)
    } else {
      Sigma <- initial
    }

    delta <- 1e-04

    Sigma.old <- matrix(0, p, p)
    count2 <- 0
    while (max(abs(Sigma - Sigma.old)) > delta) {
      # loop 1: Sigma convergence

      Sigma.old <- Sigma
      for (i in 1:p) {
        # loop 2: update each row/column of Sigma

        Omega11 <- solve(Sigma[-i, -i])
        beta <- Sigma[-i, i]

        S11 <- samp_cov[-i, -i]
        s12 <- samp_cov[-i, i]
        s22 <- samp_cov[i, i]

        a <- t(beta) %*% Omega11 %*% S11 %*% Omega11 %*% beta - 2 * t(s12) %*% Omega11 %*% beta + s22

        if (rho == 0) {
          gamma <- a
        } else if (c(a) < 10^-10) {
          gamma <- a
        } else {
          gamma <- (-1 / (2 * rho) + (1 / (4 * rho^2) + c(a) / rho)^0.5)
        }

        V <- Omega11 %*% S11 %*% Omega11 / gamma + rho * Omega11
        u <- t(s12) %*% Omega11 / gamma

        beta.old <- 0
        while (max(abs(beta - beta.old)) > delta) {
          # loop 3: off-diagonals convergence
          beta.old <- beta
          for (j in 1:(p - 1)) {
            # loop 4: each element
            temp <- u[j] - V[j, -j] %*% beta[-j]
            beta[j] <- sign(temp) * max(0, abs(temp) - rho) / V[j, j]
          }  # loop 4
        }  # loop 3

        Sigma[i, -i] <- t(beta)
        Sigma[-i, i] <- beta
        Sigma[i, i] <- gamma + t(beta) %*% Omega11 %*% beta

      }  # loop 2

      # record spcov iterations
      count2 <- count2 + 1

      if (count2 > 100) {
        cat("Omega0 fails to converge for lam1 =", rho, "lam2 =", lambda2 / K, "lam3 =", lambda3)
        break
      }

    }  # loop 1

    return(Sigma)

  }  # end of function

  # Inputs:
  K <- length(x)
  p <- dim(x[[1]])[2]

  Sa <- sapply(x, cov, simplify = "array")
  Sl <- lapply(x, cov)

  # Initial values
  Omega0 <- solve(apply(Sa, 1:2, mean) + diag(1, p) * 0.01)

  Omegas <- sapply(Sl, function(x1) solve(x1 + diag(1, p) * 0.01), simplify = "array")


  Omega0.old <- matrix(0, p, p)
  Omegas.old <- array(0, c(p, p, K))
  count <- 0

  rho <- lambda1 / (1 + lambda2 / K)

  # Start BCD algorithm
  while (max(abs(Omega0 - Omega0.old)) > delta | max(abs(Omegas - Omegas.old)) > delta) {

    # Exit if exceeds max.iters
    if (count > max.iters) {
      cat("Failed to converge for lambda1 =", lambda1, ", lambda2 =", lambda2, ", lambda3 =", lambda3,
          "delta0 =", max(abs(Omega0 - Omega0.old)),
          ", deltaK =", max(abs(Omegas - Omegas.old)))

      # Returning results
      res <- list(Omega0, Omegas)
      names(res) <- c("Omega0", "Omegas")
      return(res)
    }

    # Record current Omega0 & Omegas
    Omega0.old <- Omega0
    Omegas.old <- Omegas

    # 1st step:

    sk <- sapply(Sl, function(x1) (solve(Omega0) * lambda2 / K + x1) / (1 + lambda2 / K), simplify = FALSE)
    for (k in 1:K) {
      Omegas[, , k] <- glasso::glasso(sk[[k]], rho, penalize.diagonal = FALSE)$wi
    }

    # 2nd step:
    if (lambda3 == 0) {
      Omega0 <- apply(Omegas, 1:2, mean)
    } else {
      s0 <- apply(Omegas, 1:2, mean)
      log <- capture.output({
        Omega0 <- spcov_bcd(s0, lambda3, lambda2 = lambda2, lambda3 = lambda3)
      })
    }

    # Record BCD iterations
    count <- count + 1
  }

  res <- list(Omega0, Omegas)
  names(res) <- c("Omega0", "Omegas")
  return(res)
}
