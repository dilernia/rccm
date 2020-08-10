#' Gap Statistic
#'
#' This function selects the optimal number of clusters for the
#'  Random Covariance Clustering Model (RCCM) based on a Gap statistic
#'  as proposed by Tibshirani et al. (2001).
#' @param x List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}
#' @param gMax Maximum number of clusters or groups to consider. Must be at least 2.
#' @param B Number of reference data sets to generate.
#' @param zs \eqn{K} x \eqn{gMax} matrix with estimated cluster memberships for
#' each number of clusters considered.
#' @param optLambdas Data frame with 4 columns (lambda1, lambda2, lambda3, and \eqn{G}). The first 3 columns
#' are the tuning parameter values to implement the RCCM for a given number of clusters, and
#' the \eqn{G} column is the number of clusters that must range from 2 to \eqn{gMax}
#' @param ncores Number of computing cores to use if desired to run in parallel. Optional.
#' @return A list of length 3 containing:
#' \enumerate{
#' \item The optimally selected number of clusters (nclusts).
#' \item The \eqn{gMax} observed Gap statistics (gaps).
#' \item The \eqn{gMax-1} adjusted standard deviations of the simulated gap statistics (sigmas).
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
#' set.seed(1994)
#' myData <- rccSim(G = 2, clustSize = 10, p = 10, n = 177, overlap = 0.20, rho = 0.10)
#'
#' # Analyze simulated data with RCCM
#' optLambdas <- data.frame(lambda1 = 10, lambda2 = 50, lambda3 = 0.10, G = 2:3)
#' result2 <- rccm(x = myData$simDat, lambda1 = optLambdas$lambda1[1],
#'                 lambda2 = optLambdas$lambda2[1], lambda3 = optLambdas$lambda3[1],
#'                 nclusts = 2)
#' result3 <- rccm(x = myData$simDat, lambda1 = optLambdas$lambda1[2],
#'                 lambda2 = optLambdas$lambda2[2], lambda3 = optLambdas$lambda3[2],
#'                 nclusts = 3)
#'
#' # Estimated cluster memberships
#' zHats <- cbind(apply(result2$weights, MARGIN = 2, FUN = which.max),
#'                apply(result3$weights, MARGIN = 2, FUN = which.max))
#'
#' # Selecting number of clusters
#' clustRes <- gapSelect(x = myData$simDat, gMax = 3, B = 50, zs = zHats,
#' optLambdas = optLambdas)
#'
#' @export
#'
#' @references
#' Tibshirani, Robert, et al. "Estimating the Number of Clusters in a Data Set via the Gap
#' Statistic." Journal of the Royal Statistical Society: Series B (Statistical Methodology),
#' vol. 63, no. 2, 2001, pp. 411-423., doi:10.1111/1467-9868.00293.
gapSelect <- function(x, gMax, B = 100, zs, optLambdas, ncores = 1) {

  # Cleaning column names
  colnames(optLambdas)[1:3] <- tolower(colnames(optLambdas)[1:3])
  colnames(optLambdas)[4] <- toupper(colnames(optLambdas)[4])

  if(all.equal(colnames(optLambdas), c("lambda1", "lambda2", "lambda3", "G")) == FALSE) {
    stop("optLambdas must be a data frame with columns lambda1, lambda2, lambda3, and G.")
  }

  # Number of subjects, variables, and average number of observations
  K <- length(x)
  p <- ncol(x[[1]])
  n <- round(mean(sapply(x, FUN = nrow)))

  # Standardize data
  x <- lapply(x, FUN = scale)

  # Calculate sample precision matrices
  omegaks <- lapply(x, FUN = function(xk) {
    glasso::glasso(cov(xk), rho = 1e-16, penalize.diagonal = FALSE)$wi})

  # Calculate within-cluster variability measures for given clusterings
  Vs <- numeric(gMax)
  for(g in 2:gMax) {
    # Weights based on number of subjects in each cluster
    ws <- sapply(1:g, FUN = function(clust) {sum(zs[, g-1] == clust)})
    Vs[g] <- sum(sapply(1:g, FUN = function(clust) {
      inds <- which(zs[, g-1] == clust)
      if(length(inds) > 1) {
        return(apply(simplify2array(omegaks[inds]),
                     MARGIN = 1:2, FUN = var))
      } else {
        return(matrix(0, nrow = p, ncol = p))
      }}) %*% ws) / K
  }

  # Variability for only 1 cluster
  Vs[1] <- mean(apply(simplify2array(omegaks), MARGIN = 1:2, FUN = var))

  # Calculating range of observed values for precision matrices
  utInds <- upper.tri(diag(p), diag = T)
  q <- choose(p, 2)
  minMat <- apply(simplify2array(omegaks), MARGIN = 1:2, FUN = min)[utInds]
  maxMat <- apply(simplify2array(omegaks), MARGIN = 1:2, FUN = max)[utInds]

  # Instantiating matrix for within-cluster dispersion measures
  vMat <- matrix(NA, nrow = gMax, ncol = B)

  if (ncores > 1) {
    `%dopar%` <- foreach::`%dopar%`
    cl <- parallel::makeCluster(ncores) # creates a cluster with <ncore> cores
    doParallel::registerDoParallel(cl) # register the cluster
    vMat <- simplify2array(foreach::foreach(b = 1:B) %dopar% {
      # Generating B reference data sets and calculating within-cluster variabilities
      vMat <- rep(NA, times = gMax)
      while(sum(is.na(vMat)) > 0) {
        try(silent = TRUE, expr = {
          # Generating Data
          refDats <- lapply(1:K, FUN = function(k) {
            pMat <- matrix(0, ncol = p, nrow = p)
            while(min(eigen(pMat, symmetric = TRUE)$values) <= 0) {
              pMat <- matrix(0, ncol = p, nrow = p)
              pMat[utInds] <- sapply(1:(q + p), FUN = function(i) {
                runif(n = 1, min = minMat[i], max = maxMat[i])})
              pMat <- pMat + t(pMat)
              diag(pMat) <- diag(pMat) / 2
              eVals <- eigen(pMat, symmetric = TRUE)$values
              if(min(eVals) <= 0) {
                pMat <- pMat + diag(abs(rep(min(eVals), times = p)) + 0.01)
              }
            }
            return(mvtnorm::rmvnorm(n = n, sigma = chol2inv(chol(pMat))))})

          # Calculating MLE precision matrices
          pHats <- lapply(refDats, FUN = function(x){chol2inv(chol(cov(x)))})

          # Analyzing with varying number of clusters
          for(g in 2:gMax) {
            res <- rccm::rccm(x = refDats, lambda1 = optLambdas$lambda1[g - 1],
                             lambda2 = optLambdas$lambda2[g - 1], lambda3 = optLambdas$lambda3[g - 1],
                             nclusts = g)

            # Estimated cluster memberships
            zHats <- apply(res$weights, MARGIN = 2, FUN = which.max)

            # Weights based on number of subjects in each cluster
            ws <- sapply(1:g, FUN = function(clust) {sum(zHats == clust)})

            # Calculating within-cluster dispersion
            vMat[g] <- sum(sapply(1:g, FUN = function(clust) {
              inds <- which(zHats == clust)
              if(length(inds) > 0) {
                return(apply(simplify2array(pHats[inds]),
                             MARGIN = 1:2, FUN = var))
              } else {return(matrix(0, nrow = p, ncol = p))}}) %*% ws) / K
          }

          # Within-cluster variability specifying 1 cluster
          vMat[1] <- mean(apply(simplify2array(pHats), MARGIN = 1:2, FUN = var))

        })
      }
      return(vMat)
    })
    parallel::stopCluster(cl)
  } else {
    # Generating B reference data sets and calculating within-cluster variabilities
    for(b in 1:B) {
      while(sum(is.na(vMat[, b])) > 0) {
        try(silent = TRUE, expr = {
          # Generating Data
          refDats <- lapply(1:K, FUN = function(k) {
            pMat <- matrix(0, ncol = p, nrow = p)
            while(min(eigen(pMat, symmetric = TRUE)$values) <= 0) {

              pMat <- matrix(0, ncol = p, nrow = p)

              pMat[utInds] <- sapply(1:(q + p), FUN = function(i) {
                runif(n = 1, min = minMat[i], max = maxMat[i])})

              pMat <- pMat + t(pMat)

              diag(pMat) <- diag(pMat) / 2

              eVals <- eigen(pMat, symmetric = TRUE)$values

              if(min(eVals) <= 0) {
                pMat <- pMat + diag(abs(rep(min(eVals), times = p)) + 0.01)
              }
            }
            return(mvtnorm::rmvnorm(n = n, sigma = chol2inv(chol(pMat))))})

          # Calculating MLE precision matrices
          pHats <- lapply(refDats, FUN = function(x){chol2inv(chol(cov(x)))})

          # Analyzing with varying number of clusters
          for(g in 2:gMax) {
            res <- rccm::rccm(x = refDats, lambda1 = optLambdas$lambda1[g - 1],
                             lambda2 = optLambdas$lambda2[g - 1], lambda3 = optLambdas$lambda3[g - 1],
                             nclusts = g)

            # Estimated cluster memberships
            zHats <- apply(res$weights, MARGIN = 2, FUN = which.max)

            # Weights based on number of subjects in each cluster
            ws <- sapply(1:g, FUN = function(clust) {sum(zHats == clust)})

            # Calculating within-cluster dispersion
            vMat[g, b] <- sum(sapply(1:g, FUN = function(clust) {
              inds <- which(zHats == clust)
              if(length(inds) > 0) {
                return(apply(simplify2array(pHats[inds]),
                             MARGIN = 1:2, FUN = var))
              } else {return(matrix(0, nrow = p, ncol = p))}}) %*% ws) / K
          }

          # Within-cluster variability specifying 1 cluster
          vMat[1, b] <- mean(apply(simplify2array(pHats), MARGIN = 1:2, FUN = var))

        })
      }
    }
  }

  # Calculating gap statistics
  vBars <- apply(log(vMat), MARGIN = 1, FUN = mean)
  sigmas <- apply(log(vMat), MARGIN = 1, FUN = sd)*sqrt((n-1) / n)*sqrt(1 + 1/B)
  gaps <- vBars - log(Vs)

  # Selecting optimal number of clusters
  checks <- c(sapply(2:(gMax-1), FUN = function(g) {
    gaps[g-1] >= gaps[g] - sigmas[g]}), TRUE)


  # Selecting optimal number of clusters
  checks <- c(sapply(1:(gMax-1), FUN = function(g) {
    gaps[g] >= gaps[g + 1] - sigmas[g + 1]}), TRUE)

  nclusts <- (1:gMax)[which(checks)[1]]

  # Returning optimal number of clusters
  return(list(nclusts = nclusts,
              gaps = gaps, sigmas = sigmas))
}
