#' Matrix Clustering
#'
#' This function provides the probability density function for
#' the Wishart distribution.
#' @param mats List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}.
#' @param G Number of groups or clusters.
#' @return \eqn{G} x \eqn{K} matrix of weights
#'
MatClust <- function(mats, G) {
  K <- dim(mats)[3]
  combos <- expand.grid(s1 = 1:K, s2 = 1:K)
  distMat <- matrix(NA, nrow = K, ncol = K)
  for(r in 1:K) {
    for(s in 1:K) {
      distMat[r, s] <- norm(mats[, , r] - mats[, , s], type = 'F')
    }
  }

  cl0 <- cutree(hclust(d = as.dist(distMat), method = "ward.D"), k = G)

  wgk <- matrix(NA, nrow = G, ncol = K)

  for(i in 1:G) {
    for(j in 1:K) {
      wgk[i, j] <- ifelse(cl0[j] == i, 1, 0)
    }
  }
  return(wgk)
}

#' Negative Log-Likelihood
#'
#' This function calculates the negative log-likelihood for
#' selecting the optimal tuning parameters using cross-validation.
#' @param omegaks \eqn{p} x \eqn{p} x \eqn{K} array of \eqn{K} number of estimated subject-level precision matrices.
#' @param datf List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}.
#' @return Negative log-likelihood
#'
nllikCalc <- function(omegaks, datf) {
  K <- length(datf)
  p <- ncol(omegaks[, , 1])
  nlogLikes <- sum(sapply(X = 1:K, FUN = function(subj) {
    nk <- nrow(datf[[subj]])
    nllik <- sum(diag(nk*cov(datf[[subj]]) %*% omegaks[, , subj])) - nk*log(det(omegaks[, , subj]))
    return(nllik)}))
  return(nlogLikes)
}

#' Tuning Parameter Selection via Cross-Validation
#'
#' This function selects the optimal set of tuning parameters
#'  for the Random Covariance Clustering Model (RCCM) based on
#'  cross-validation.
#' @param x List of \eqn{K} data matrices each of dimension \eqn{n_k} x \eqn{p}
#' @param G Number of groups or clusters.
#' @param lambs A data frame of candidate tuning parameter values with three columns: lambda1, lambda2, and lambda3.
#' @param methods Methods to implement cross-validation for. Must be a vector containing one or more of "FGL", "GGL", "GLasso", "RCCM", or "RCM".
#' @param folds Number of folds to use for cross-validation.
#'
#' @return A data frame of optimally selected tuning parameter values and the
#' correspond average negative log-likelihoods across the \eqn{folds} number
#' of folds with five columns: lambda1, lambda2, lambda3, nll, and method.
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
#' G <- 2
#' myData <- rccSim(G = G, clustSize = 10, p = 10, n = 177, overlap = 0.20, rho = 0.10)
#'
#' # Analyze simulated data with RCCM & GLasso using 5-fold CV
#' lambdas <- expand.grid(lambda1 = c(5, 10, 15), lambda2 = c(50, 100), lambda3 = 0.10)
#' cvLambdas <- cvTune(x = myData$simDat, G = G, lambs = lambdas, methods = c("RCCM", "GLasso"), folds = 5)
#'
#' # RCCM with CV selected tuning parameters
#' resRCCM <- rccm(x = myData$simDat, lambda1 = cvLambdas[which(cvLambdas$method == "RCCM"), ]$lambda1,
#'                 lambda2 = cvLambdas[which(cvLambdas$method == "RCCM"), ]$lambda2,
#'                 lambda3 = cvLambdas[which(cvLambdas$method == "RCCM"), ]$lambda3, nclusts = 2)
#'
#' # GLasso with CV selected tuning parameters
#' resGLasso <- lapply(myData$simDat, FUN = function(x) {
#' glasso::glasso(cov(x), rho = cvLambdas[which(cvLambdas$method == "GLasso"), ]$lambda1,
#' penalize.diagonal = FALSE)$wi})
#'
#' @export
cvTune <- function(x, G, lambs, methods = c("RCCM", "GGL", "GLasso", "FGL", "RCM"), folds = 5) {

  # Cleaning column names
  colnames(lambs)[1:3] <- tolower(colnames(lambs)[1:3])

  if(all.equal(colnames(lambs), c("lambda1", "lambda2", "lambda3")) == FALSE) {
    stop("lambs must be a data frame with columns lambda1, lambda2, and lambda3")
  }

  # Number of subjects and tuning parameter combinations
  K <- length(x)
  J <- nrow(lambs)

  # Randomly reorder rows of data
  simDatsCV <- lapply(x, FUN = function(datf) {return(datf[sample(nrow(datf)), ])})

  # Storing indices for each fold
  keeps <- cut(seq(1, nrow(simDatsCV[[1]])), breaks = folds, labels = FALSE)

  # Calculating negative log likelihood for each fold for each tuning parameter set
  nll <- t(sapply(1:J, simplify = "array", FUN = function(j) {
    apply(sapply(X = 1:folds, FUN = function(fold, datf = simDatsCV, lambdas = lambs) {

      # Instantiating list of results
      resList <- list()

      # Separating into training and test sets
      xTest <- lapply(datf, FUN = function(x) {return(x[which(keeps == fold), ])})
      datf <- lapply(datf, FUN = function(x) {return(x[which(keeps != fold), ])})

      # Standardizing test data
      xTest <- lapply(xTest, FUN = scale)

      if("RCCM" %in% methods) {
        # Conducting analysis of simulated data using RCCM
        RCCMres <- list()
        t1 <- Sys.time()
        RCCMres$results <- rccm(x = datf, lambda1 = lambdas$lambda1[j],
                                lambda2 = lambdas$lambda2[j], lambda3 = lambdas$lambda3[j],
                                nclusts = G)
        RCCMres$Time <- Sys.time() - t1

        # Adding to results list
        resList[[length(resList)+1]] <- RCCMres
        names(resList)[length(resList)] <- "RCCM"
      }

      K <- length(datf)
      Sl <- sapply(datf, cov, simplify = "array")
      gl <- sapply(1:K, simplify = "array",
                   FUN = function(x){glasso::glasso(Sl[, , x], rho = 0.001,
                                                    penalize.diagonal = FALSE)$wi})

      if("FGL" %in% methods) {
        # Conducting analysis using Hierarchical Clustering then Fused Lasso
        FGLres <- list()

        # Initializing weights using hierarchical clustering based on dissimilarity matrix of Frob norm differences
        FGLres$results$weights <- MatClust(gl, G = G)

        # Fused lasso penalty
        t1 <- Sys.time()
        jglList <- unlist(lapply(FUN = function(g) {
          if(length(which(as.logical(FGLres$results$weights[g, ]))) > 1) {
            prec <- JGL::JGL(Y = datf[which(as.logical(FGLres$results$weights[g, ]))], penalty = "fused",
                             penalize.diagonal=FALSE, lambda1 = lambdas$lambda1[j]/100,
                             lambda2 = lambdas$lambda2[j]/50/1000,
                             return.whole.theta = TRUE)$theta
          } else {
            prec <- list(glasso::glasso(s = cov(datf[which(as.logical(FGLres$results$weights[g, ]))][[1]]),
                                        rho = lambdas$lambda1[j]/100, penalize.diagonal = FALSE)$wi)
          }
          return(setNames(prec, c(which(as.logical(FGLres$results$weights[g, ])))))
        }, X = 1:G), recursive = F)
        FGLres$results$Omegas <- jglList[as.character(1:length(jglList))]
        FGLres$Time <- Sys.time() - t1

        # Adding to results list
        resList[[length(resList)+1]] <- FGLres
        names(resList)[length(resList)] <- "FGL"
      }

      if("GGL" %in% methods) {
        # Conducting analysis using Hierarchical Clustering then Group Lasso
        GGLres <- list()

        # Initializing weights using hierarchical clustering based on dissimilarity matrix of Frob norm differences
        GGLres$results$weights <- MatClust(gl, G = G)

        # Group lasso penalty
        t1 <- Sys.time()
        jglList <- unlist(lapply(FUN = function(g) {
          if(length(which(as.logical(GGLres$results$weights[g, ]))) > 1) {
            prec <- JGL::JGL(Y = datf[which(as.logical(GGLres$results$weights[g, ]))], penalty = "group",
                             penalize.diagonal=FALSE, lambda1 = lambdas$lambda1[j]/100,
                             lambda2 = lambdas$lambda2[j]/50/1000, return.whole.theta = TRUE)$theta
          } else {
            prec <- list(glasso::glasso(s = cov(datf[which(as.logical(GGLres$results$weights[g, ]))][[1]]),
                                        rho = lambdas$lambda1[j]/100, penalize.diagonal = FALSE)$wi)
          }
          return(setNames(prec, c(which(as.logical(GGLres$results$weights[g, ])))))
        }, X = 1:G), recursive = F)
        GGLres$results$Omegas <- jglList[as.character(1:length(jglList))]
        GGLres$Time <- Sys.time() - t1

        # Adding to results list
        resList[[length(resList)+1]] <- GGLres
        names(resList)[length(resList)] <- "GGL"
      }

      if("GLasso" %in% methods) {
        # Conducting analysis using Hierarchical Clustering then Group Lasso
        GLassores <- list()

        # Initializing weights using hierarchical clustering based on dissimilarity matrix of Frob norm differences
        GLassores$results$weights <- MatClust(gl, G = G)

        # Individual GLasso after hierarchical clustering
        t1 <- Sys.time()
        GLassores$results$Omegas <- lapply(FUN = function(x) {
          return(glasso::glasso(Sl[, , x], rho = lambdas$lambda1[j]/100, penalize.diagonal = FALSE)$wi)
        }, X = 1:(dim(Sl)[3]))
        GLassores$Time <- Sys.time() - t1

        # Adding to results list
        resList[[length(resList)+1]] <- GLassores
        names(resList)[length(resList)] <- "GLasso"
      }

      if("RCM" %in% methods) {
        # Conducting analysis using Hierarchical Clustering then Group Lasso
        RCMres <- list()

        # Initializing weights using hierarchical clustering based on dissimilarity matrix of Frob norm differences
        RCMres$results$weights <- MatClust(gl, G = G)

        # Individual GLasso after hierarchical clustering
        t1 <- Sys.time()
        rcmList <- lapply(FUN = function(g) {
          rcmTemp <- randCov(x = datf[which(as.logical(RCMres$results$weights[g, ]))],
                             lambda1 = lambdas$lambda1[j]/100,
                             lambda2 = lambdas$lambda2[j]/50,
                             lambda3 = lambdas$lambda3[j]/100000)
          prec <- rcmTemp$Omegas
          Omega0 <- rcmTemp$Omega0
          return(list(Omegas = setNames(prec, c(which(as.logical(RCMres$results$weights[g, ])))),
                      Omega0 = Omega0))
        }, X = 1:G)

        # Storing arrays for subject and group-level estimates
        RCMres$results$Omegas <- sapply(1:K, FUN = function(x){
          for(g in 1:length(rcmList)) {
            if(as.character(x) %in% names(rcmList[[g]]$Omegas)) {
              return(rcmList[[g]]$Omegas[, ,  which(names(rcmList[[g]]$Omegas) == as.character(x))])
            }}}, simplify = "array")

        RCMres$results$Omega0 <- sapply(1:G, FUN = function(x){
          return(rcmList[[x]]$Omega0)}, simplify = "array")

        RCMres$Time <- Sys.time() - t1

        # Adding to results list
        resList[[length(resList)+1]] <- RCMres
        names(resList)[length(resList)] <- "RCM"
      }

      # Creating vector of negative log-likelihoods from test data
      nlogLikes <- c()
      for(m in methods){
        if(is.list(resList[[m]]$results$Omegas)) {
          Oms <- sapply(resList[[m]]$results$Omegas, FUN = function(omegas){return(omegas)}, simplify = "array")
        } else {
          Oms <- resList[[m]]$results$Omegas
        }
        nlogLikes <- c(nlogLikes, nllikCalc(omegaks = Oms, datf = xTest))
      }

      return(nlogLikes)
    }),
    FUN = mean, MARGIN = 1)}))

  nmethods <- length(methods)

  cvRes <- data.frame(lambda1 = rep(lambs$lambda1, times = nmethods),
                      lambda2 = rep(lambs$lambda2, times = nmethods),
                      lambda3 = rep(lambs$lambda3, times = nmethods),
                      nll = as.numeric(nll), method = rep(methods, each = J))

  return(do.call(rbind, by(cvRes, cvRes$method, function(x) {x[which.min(x$nll), ]})))
}
