#' Modified Stability Approach for Regularization Selection
#'
#' This function implements a modified stability approach for
#' regularization selection (stARS) method for tuning parameter
#' selection. Methods available to implement include the
#' fused graphical lasso (\link[JGL:JGL]{FGL}), group graphical lasso (\link[JGL:JGL]{GGL}),
#' graphical lasso (\link[glasso:glasso]{GLasso}), random covariance clustering
#' model (\link[rcm:rccm]{RCCM}), and the random covariance model (\link[rcm:randCov]{RCM}).
#'
#' @param datf List of \eqn{K} data sets each of dimension \eqn{n_k} x \eqn{p}.
#' @param lambs A data frame of candidate tuning parameter values with three columns: lambda1, lambda2, and lambda3.
#' @param method Method to implement modified stARS algorithm for. Must be one of "FGL", "GGL", "GLasso", "RCCM", or "RCM".
#' @param G Number of groups or clusters. Only applicable if method = "RCCM".
#' @param N Number of subsamples for modified stARS algorithm
#' @param beta Positive scalar between 0 and 1. Limits allowed amount of instability across subsamples.
#' @param z0s Vector of length \eqn{K} with initial cluster memberships. Only applicable if method = "RCCM".
#' @param ncores Number of computing cores to use if desired to run in parallel. Optional.
#' @return A data frame of optimally selected tuning parameter values and the sparsity level with four columns: lambda1, lambda2, lambda3, and sparsity.
#'
#' @author
#' Andrew DiLernia
#'
#' @examples
#' # Generate data with 2 clusters with 10 subjects in each group,
#' # 10 variables for each subject, 100 observations for each variable for each subject,
#' # the groups sharing about 50% of network connections, and 10% of differential connections
#' # within each group
#' set.seed(1994)
#' myData <- rccSim(G = 2, clustSize = 10, p = 10, n = 100, overlap = 0.50, rho = 0.10)
#'
#' # Find optimal tuning parameter set using modified stARS
#' optTune <- starsRccm(datf = myData$simDat, lambs = expand.grid(lambda1 = c(20, 25, 30),
#' lambda2 = c(300, 325), lambda3 = 0.01), method = "RCCM", G = 2)
#'
#' # Analyze with RCCM using optimally selected tuning parameters
#' resultRccm <- rccm(x = myData$simDat, lambda1 = optTune$lambda1[1],
#' lambda2 = optTune$lambda2[1], lambda3 = optTune$lambda3[1], nclusts = 2)
#'
#' @export
#'
#' @references
#' Liu, Han, Kathryn Roeder, and Larry Wasserman.
#' "Stability Approach to Regularization Selection (StARS) for High Dimensional Graphical Models." 2010.
#'
#' Danaher, Patrick, Pei Wang, and Daniela M. Witten.
#' "The Joint Graphical Lasso for Inverse Covariance Estimation across Multiple Classes." Journal of the Royal Statistical Society: Series B (Statistical Methodology) 76, no. 2 (2014): 373-97.
#'
#' Friedman, Jerome, Trevor Hastie, and Robert Tibshirani.
#' "Sparse Inverse Covariance Estimation with the Graphical Lasso." Biostatistics 9, no. 3 (2008): 432-41.
#'
#' Zhang, Lin, Andrew DiLernia, Karina Quevedo, Jazmin Camchong, Kelvin Lim, and Wei Pan.
#' "A Random Covariance Model for Bi-level Graphical Modeling with Application to Resting-state FMRI Data." 2019. https://arxiv.org/pdf/1910.00103.pdf
starsRccm <- function(datf, lambs, method = "RCCM", G = 2, N = 10,
                      beta = 0.05, z0s = NULL, ncores = NULL) {

  ns <- sapply(datf, FUN = nrow)
  K <- length(datf)

  # Setting b: size of subsamples, N: number of subsamples to draw, beta: ceiling for instability measure
  bs <- sapply(ns, FUN = function(x) {
    if (x > 100) {
      floor(10 * sqrt(x))
      }
    else {
      floor(0.75 * (x))}})

  # Drawing subsamples each of size b from data until
  # N unique subsamples have been collected
  keeper <- function(b, n) {
    keepInds <- matrix(NA, nrow = b, ncol = N)
    while (ncol(unique(keepInds, MARGIN = 2)) < N) {
      for (s in 1:N) {
        keepInds[, s] <- sort(sample(x = 1:n, size = b, replace = FALSE))
      }
    }
    return(keepInds)
  }

  keepInds <- mapply(FUN = keeper, n = ns, b = bs, SIMPLIFY = FALSE)

  # Reducing lambdas when possible
  if (method == "GLasso") {
    lambs <- data.frame(lambda1 = unique(lambs$lambda1), lambda2 = NA, lambda3 = NA)
  } else if (method %in% c("FGL", "GGL")) {
    lambs <- expand.grid(lambda1 = unique(lambs$lambda1), lambda2 = unique(lambs$lambda2), lambda3 = NA)
  } else if (method %in% c("RCCM", "RCM")) {
    lambs <- lambs
  }

  # Updating number of needed cores
  if(is.null(ncores)) {
    ncores <- 1
  }
  ncores <- min(c(ncores, nrow(lambs)))

  # Implement select method for each of N bootstrap subsamples and obtaining networks
  starNets <- lapply(1:N, FUN = function(i) {

    # Obtaining ith bootstrap sample of data
    subDats <- lapply(1:K, FUN = function(k) {
      datf[[k]][keepInds[[k]][, i], ]})

    # Running rccm method for bootstrap sample for each lambda combination in parallel if requested
    if (ncores > 1) {
      `%dopar%` <- foreach::`%dopar%`
      cl <- parallel::makeCluster(ncores) # creates a cluster with <ncore> cores
      doParallel::registerDoParallel(cl) # register the cluster
      nets <- foreach::foreach(t = 1:nrow(lambs),
                               .export = c("method", "G", "K", "lambs", "z0s")) %dopar% {
                                 listRes <- NULL
                                 tryCatch({
                                   if (method == "RCCM") {
                                     arrayRes <- rccm(subDats, lambda1 = lambs[t, "lambda1"], lambda2 = lambs[t, "lambda2"],
                                                      lambda3 = lambs[t, "lambda3"], nclusts = G, z0s = z0s)$Omegas
                                     listRes <- lapply(lapply(1:K, FUN = function(k) {
                                       arrayRes[, , k]}), FUN = adj)
                                                                        } else if (method == "GLasso") {
                                     listRes <- lapply(subDats, FUN = function(x) {
                                       adj(glasso::glasso(cov(x), rho = lambs[t, "lambda1"] / 100, penalize.diagonal = FALSE)$wi)})
                                                                        } else if (method %in% c("GGL", "FGL")) {
                                     listRes <- lapply(JGL::JGL(Y = subDats, penalty = ifelse(method == "GGL", "group", "fused"),
                                                                penalize.diagonal = FALSE,
                                                                lambda1 = lambs[t, "lambda1"] / 100,
                                                                lambda2 = lambs[t, "lambda2"] / 50 / 1000,
                                                                return.whole.theta = TRUE)$theta, FUN = adj)
                                                                        } else if (method == "RCM") {
                                     arrayRes <- randCov(x = subDats, lambda1 = lambs[t, "lambda1"] / 100,
                                                                    lambda2 = lambs[t, "lambda2"] / 50,
                                                                    lambda3 = lambs[t, "lambda3"] / 100000)$Omegas
                                     listRes <- lapply(lapply(1:K, FUN = function(k) {
                                       arrayRes[, , k]}), FUN = adj)
                                                                        }
                                 }, error = function(e) {
                                   warning(paste0("stARS failed for lambda1 = ", lambs[t, "lambda1"], ", lambda2 = ",
                                                  lambs[t, "lambda2"], ", lambda3 = ", lambs[t, "lambda3"]))
                                   return(NULL)
                                 })
                                 return(listRes)
                               }
      parallel::stopCluster(cl)
    } else {
      nets <- lapply(1:nrow(lambs), FUN = function(t) {
        tryCatch({
          if (method == "RCCM") {
            arrayRes <- rccm::rccm(subDats, lambda1 = lambs[t, "lambda1"], lambda2 = lambs[t, "lambda2"],
                             lambda3 = lambs[t, "lambda3"], nclusts = G, z0s = z0s)$Omegas
            listRes <- lapply(lapply(1:K, FUN = function(k) {
              arrayRes[, , k]}), FUN = adj)
          } else if (method == "GLasso") {
            listRes <- lapply(subDats, FUN = function(x) {
              adj(glasso::glasso(cov(x), rho = lambs[t, "lambda1"] / 100,
                                                                            penalize.diagonal = FALSE)$wi)})
          } else if (method %in% c("GGL", "FGL")) {
            listRes <- lapply(JGL::JGL(Y = subDats, penalty = ifelse(method == "GGL", "group", "fused"),
                                       penalize.diagonal = FALSE,
                                       lambda1 = lambs[t, "lambda1"] / 100,
                                       lambda2 = lambs[t, "lambda2"] / 50 / 1000,
                                       return.whole.theta = TRUE)$theta, FUN = adj)
          } else if (method == "RCM") {
            arrayRes <- randCov(x = subDats, lambda1 = lambs[t, "lambda1"] / 100,
                                           lambda2 = lambs[t, "lambda2"] / 50,
                                           lambda3 = lambs[t, "lambda3"] / 100000)$Omegas
            listRes <- lapply(lapply(1:K, FUN = function(k) {
              arrayRes[, , k]}), FUN = adj)
          }
          return(listRes)
        }, error = function(e) {
          warning(paste0("stARS failed for lambda1 = ", lambs[t, "lambda1"], ", lambda2 = ",
                         lambs[t, "lambda2"], ", lambda3 = ", lambs[t, "lambda3"]))
          return(NULL)
          })
      })
    }
    return(nets)
  })

  # Function for calculating total instability statistic (D stat) for each lambda
  dCalc <- function(lambda, aMats = starNets) {

    # Subsetting to only include matrices for input lambda value
    aMats <- lapply(aMats, FUN = function(m) {
      m[[lambda]]})

    # Calculating theta and E matrices
    thetaMats <- lapply(1:K, FUN = function(k) {
      Reduce("+", lapply(aMats, FUN = function(m) {
        m[[k]]})) /
        sum(sapply(aMats, FUN = function(m) {
          !is.null(m[[k]])}))})

    eStatMats <- lapply(thetaMats, FUN = function(m) {
      2 * m * (1 - m)})

    # Calculating total instability statistic (D stat)
    Dstat <- mean(sapply(eStatMats, FUN = function(x) {
      if (!is.null(dim(x))) {
        Dstat <- 0

        for (r in 1:nrow(x)) {
          for (c in 1:ncol(x)) {
            if (r < c) {
              Dstat <- Dstat + x[r, c]
            }
          }
        }
        return(Dstat)
      } else {
        NULL
        }
    }))

    # Calculating sparsity level
    sparsity <- 1 - mean(sapply(thetaMats, FUN = function(m) {
      unlist(m[lower.tri(m, diag = FALSE)])}))

    if (!is.null(ncol(eStatMats[[1]]))) {
      results <- data.frame(D = Dstat / choose(n = ncol(eStatMats[[1]]), k = 2),
                            Sparsity = sparsity, Lambda = lambda)
      return(results)
    } else {
      return(data.frame(D = NA, Sparsity = NA, Lambda = lambda))
      }
  }

  dResults <- do.call("rbind", lapply(X = 1:nrow(lambs), FUN = dCalc))

  dResults$dBar <- NA

  # Calculating dBar value for each lambda
  for (r in 1:nrow(dResults)) {
    if (!is.na(dResults[r, c("D")])) {
      thresh <- as.numeric(dResults[r, "Sparsity"])
      dResults[r, c("dBar")] <- max(unlist(dResults[which(dResults$Sparsity >= thresh), c("D")]))
    }
  }

  # Calculating optimal lambda value
  if (length(which(dResults$dBar <= beta)) > 0) {
    optSparse <- min(unlist(dResults[which(dResults$dBar <= beta), c("Sparsity")]))
    optD <- min(unlist(dResults[which(dResults$dBar <= beta & dResults$Sparsity == optSparse), c("D")]))
  } else {
    # Displaying a warning that instability was never low enough
    warning(paste0("Minimum instability of ", round(min(dResults$dBar), 3), " not less than beta = ", beta,
                   ". Returning maximum sparsity result."))
    optSparse <- max(dResults$Sparsity)
    optD <- min(dResults[which(dResults$Sparsity == optSparse), "D"])
  }
  starLamb <- cbind(lambs[dResults[which(dResults$Sparsity == optSparse & dResults$D == optD),
                                   "Lambda"], ], sparsity = optSparse)

  return(starLamb)
}
