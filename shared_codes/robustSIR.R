# --------------- Description ------------------------- #
# * This is the main file containing the functions to implement ROSE algorithm for solving the penalized optimization problem (10):
#   - robustSIR: The main code of ROSE algorithm
#   - cv.robustSIR: Cross-validation function of ROSE algorithm.
#   - msda, cv.msda: The revised 'msda' and 'cv.msda' functions from package 'msda'. These two functions are used for automatically generating tuning parameter sequences in 'robustSIR' and 'cv.robustSIR' functions.
#   - MU: Calculates the Sigma and U matrices in eq (9).
#   - admm: The main code of ADMM algorithm, whose detailed description is presented in Section S.5.
#   - updataB: Update G iterate in the ADMM algorithm for solving eq. (S.5).
#   - updateC: Update C iterate in the ADMM algorithm for solving eq. (S.5).
# ---------------------------------------- #

# ---------- Main program of ROSE algorithm ---------- #
# * Inputs:
# ========
# x: n x p observation matrix for predictor.
# y: n-dimensional observation vector for response.
# yclass: Discretized response taking values in 1,...,H.
# d: True structural dimension. The default is NULL.
# categorical: A logical value indicating whether y is categorical.
# H: The number of slices. The default value is 5.
# M: The Sigma matrix in optimization problem in eq (10). If both arguments 'M' and 'U' are provided, the arguments "x" and "y" are ignored.
# U: The U matrix in optimization problem in eq (10). If both arguments 'M' and 'U' are provided, the arguments "x" and "y" are ignored.
# nobs: The number of observations.
# nlambda: the length of the tuning parameter sequence in the cv.msda function. This argument is passed to the "cv.msda" function.
# lambda.factor: the factor used for generating the tuning parameter sequence in the cv.msda function. A smaller factor gives a sparsely spaced sequence while a larger factor gives a dense sequence. This argument is passed to the "cv.msda" function.
# lam1, lam2, gamma: The user-specified sequences of tuning parameters lambda_1, lambda_2 and gamma. If NULL the sequences are automatically generated. 
# lam1_fac, lam2_fac: The factors used in automatically generating the tuning parameter sequences.
# eps: The tolerance of convergence in ADMM algorithm. This argument is passed to the "admm" function.
# maxit: The maximal iterations in ADMM algorithm. This argument is passed to the "admm" function.
# warning.it: The warning level. The higher this value is, the more warning messages will be printed. If warning.it = 0, all warning messages will be suppressed.
# true.beta: The true value of the basis matrix of the central subspace. This argument is only used in the simulation studies for checking the discrepancy between the estimated and true basis matrices.
# ...: Other arguments passed to the "admm" function
# 
# * Outputs:
# ========
# beta: A list containing the estimated basis matrices of central subspace.
# B: A list containing estimated B matrices.
# rank: A vector containing the estimated ranks.
# s: A vector containing the estimated sparsity levels.
# dist: A vector containing the subspace distances between each estimated subspace in the solution path and the true subspace, if "true.beta" is not null.
# lam1, lam2, gamma: The tuning parameter sequences.
# step: A vector containing the number of iterations to converge for each tuning parameter.
# time: The execution time of implementing the ADMM algorithm for solving eq (10) once.
# code: Error code. If code == 0, exit normally.

robustSIR <- function(x = NULL, y = NULL, yclass = NULL, d = NULL, categorical=FALSE, H=5, M = NULL, U = NULL, nobs = NULL, nlambda=NULL, lambda.factor=NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.3, length.out = 10), lam2_fac=seq(0.01,0.3, length.out = 10), eps = 1e-3, maxit = 1e+3, warning.it = 1, true.beta = NULL, ...){

    if(is.null(M) || is.null(U)){
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)){
      if(categorical == FALSE){
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
        nclass <- as.integer(length(unique(yclass)))
      }
      else if(categorical == TRUE){
        yclass <- y
      }
    }
    if(is.null(gamma)){
      gamma <- c(10,30,50)
    }
    if(is.null(lam1) || is.null(lam2)){
      fit_1 <- cv.msda(x, y, yclass = yclass, nlambda=nlambda, lambda.factor=lambda.factor, nfolds = 5, maxit=1e3, warning.it = warning.it, categorical = categorical)
      M <- fit_1$M
      U <- fit_1$U
      M_fold <- fit_1$M_fold
      U_fold <- fit_1$U_fold
      id_min_msda <- fit_1$id
      lam1_min_msda <- fit_1$lam_min
      beta_msda <- as.matrix(fit_1$beta)
      if(is.null(lam1)) lam1 <- (lam1_min_msda)*lam1_fac
      if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
      if (all(lam2 == 0)){
        lam2 <- 0
        if(warning.it >= 1) warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.")
      }
    }else{
      MU_out <- MU(x, y, yclass)
      M <- MU_out$M
      U <- MU_out$U
    }
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
  }
  else{
    if(is.null(lam1) || is.null(lam2) || is.null(gamma)) stop("Sequences lam1, lam2 or gamma is missing.")
    if(is.null(nobs)) stop("Missing nobs.")
    nvars <- NCOL(M)
  }
  
  ## Error code
  code <- 0
  
  if(is.vector(lam1) && (length(lam1) == 1) && (lam1 == 0) && is.vector(lam2) && (length(lam2) == 1) && (lam2 == 0)){
    B <- solve(M) %*% U
    if(is.null(d)) beta <- svd(B)$u
    else if(d == 0) beta <- matrix(0, nrow(Bnew), ncol(Bnew))
    else beta <- svd(B)$u[,1:d,drop=FALSE]
    vec <- as.vector(beta)
    vec[abs(vec) < 1e-3] <- 0
    beta <- matrix(vec, nrow(beta), ncol(beta))
    rank <- NCOL(beta)
    output <- list(beta = beta, B = B, rank = rank, lam1 = lam1, lam2 = lam2, code = code)
  }
  else{
    # Fit with admm function
    fit <- admm(M, U, nobs, nvars, lam1, lam2, gamma, eps, maxit, d, warning.it = warning.it, ...)
    B_l <- fit$B
    beta_l <- fit$beta
    dist_l <- rep(NA, length(beta_l))
    if(!is.null(true.beta)){
      for (idx in seq_along(beta_l)){
        if(!is.null(beta_l[[idx]])){dist_l[idx] <- subspace(beta, beta_l[[idx]])}        
      }
    }
    if (all(sapply(beta_l, is.null))){
      code <- 1
      rank_l <- rep(NA, length(beta_l))
      s_l <- rep(NA, length(beta_l))
      dist_l <- rep(NA, length(beta_l))
      if(warning.it >= 2) warning("In \"Robust SIR\": No converged results returned.\n")
      if(length(B_l) == 1){
        B_l <- B_l[[1]]
        beta_l <- beta_l[[1]]
        rank_l <- NA
        s_l <- NA
      }
      return(list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, dist = dist_l, step = NA, time = NA, code = code))
    }
    rank_l <- fit$rank
    s_l <- fit$s
    step_l <- fit$step
    time_l <- fit$time
    if(length(B_l) == 1){
      B_l = B_l[[1]]; beta_l = beta_l[[1]]; rank_l = rank_l[[1]]; s_l = s_l[[1]]; step_l = step_l[[1]]; time_l = time_l[[1]]; dist_l = dist_l[[1]]
    }
    output <- list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, dist = dist_l, lam1 = lam1, lam2 = lam2, gamma = gamma, step = step_l, time = time_l, code = code)
  }
  output
}

# ---------- Cross-validation function of ROSE algorithm ---------- #
# * The inputs and outputs are similar to those in "robustSIR" function. Only the distinct function are described in the following. Please refer to the description of the "robustSIR" function for other input/output arguments.
#
# * Inputs:
# =======
# nfolds: The number of folds in the cross-validation.
# fold: User-defined fold index.
# lambda.type: The type for selecting tuning parameters.
#   - min: selects the tuning parameters with the minimal cv-error.
#   - 1se.rank: among the tuning parameters yielding the cv-error within one-standard-error of the minimal cv-error, selects the tuning parameters corresponding to the basis matrix with the smallest rank.
#   - 1se.s: among the tuning parameters yielding the cv-error within one-standard-error of the minimal cv-error, selects the tuning parameters corresponding to the basis matrix with the smallest sparsity level.
# trace.it: If TRUE, print the process of cross-validation.
# original_y: A logical value. If TRUE, use the original y in computing the distance correlation, which is used as the tuning parameter selection criterion. Otherwise, use the discretized y in computing the distance correlation. This argument is passed to the "eval_dc" function.
# solution.path: A logical value. If TRUE, and "plot=TRUE", generates a solution path of the cv-error.
# solution.path.distance: A logical value. If TRUE, "true.beta" is not null, and "plot=TRUE", generates a solution path of the subspace distances between each estimated subspace and the true subspace.
# solution.path.s: A logical value. If TRUE, and "plot=TRUE", generates a solution path of the estimated sparsity level.
# plot: If TRUE, (1) plot a cv-error curve in the "cv.msda" function; (2) plot a cv-error curve in each fold; (3) plot a complete solution path of the cv-error if "solution.path=TRUE"; (4) plot a solution path of the subspace distances if "solution.path.distance=TRUE" and "true.beta" is not null; (5) plot a solution path of the estimated sparsity level if "solution.path.s=TRUE".
# 
# * Outputs:
# ========
# eval: Cv-errors in all folds.
# cvm: Mean of cv-errors over folds.
# id_lam1, id_lam2, id_gamma: Index of the optimal tuning parameters in the sequence.
# lam1_hat, lam2_hat, gamma_hat: The optimal tuning parameters.
# fold: The user-defined fold index.
# beta_path, rank_path, s_path: the solution paths of beta matrix, rank, and sparsity level.
# gammahat: Estimated gamma in eq (5).
# what: Estimate W in eq (7).

cv.robustSIR <- function(x, y, yclass = NULL, d = NULL, categorical=FALSE, H=5, lambda.factor=NULL, nlambda=NULL, nfolds = 5, fold = NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.3, length.out = 10), lam2_fac=seq(0.01,0.3, length.out = 10), lambda.type = c("min", "1se.rank", "1se.s"), plot = FALSE, eps = 1e-3, maxit = 1e+3, trace.it = TRUE, solution.path = FALSE, true.beta = NULL, warning.it = 1, solution.path.distance = FALSE, solution.path.s = FALSE, original_y = FALSE, ...){

  start_time <- Sys.time() # Record time
  
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  if(any(table(yclass) < nfolds)){
    stop(sprintf("The sample size of class %d is less than nfolds\n", which(table(yclass) < 5)))
  }
  
  lambda.type <- match.arg(lambda.type)
  ord <- order(y)
  y <- y[ord]
  yclass <- yclass[ord]
  x <- x[ord,]
  M <- U <- M_fold <- U_fold <- NULL
  
  nobs <- dim(x)[1]
  if (nfolds < 3) stop("nfolds must be larger than 3")
  if (nfolds > nobs) stop("nfolds is larger than the sample size")
  count <- as.numeric(table(yclass))
  
  if(is.null(fold)){
    fold <- c()
    for(cnt in count){
      fold <- c(fold, sample(rep(seq(nfolds), length = cnt)))
    }
  }
  else{
    nfolds <- length(unique(fold))
  }
  
  if(is.null(gamma)){
    gamma <- c(10,30,50)
  }
  if(is.null(lam1) || is.null(lam2)){
    fit_1 <- cv.msda(x, y, yclass = yclass, nlambda=nlambda, lambda.factor=lambda.factor, fold = fold, maxit=1e3, warning.it =  warning.it, plot = plot, categorical = categorical, original_y = original_y)
    M <- fit_1$M
    U <- fit_1$U
    what <- fit_1$what
    M_fold <- fit_1$M_fold
    U_fold <- fit_1$U_fold
    gammahat <- fit_1$gammahat
    id_min_msda <- fit_1$id
    lam1_min_msda <- fit_1$lam_min
    beta_msda <- as.matrix(fit_1$beta)
    if(is.null(lam1)) lam1 <- (lam1_min_msda)*lam1_fac
    if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
    if (all(lam2 == 0)){
      lam2 <- 0
      if(warning.it >= 1) warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.")
    }
  }
  n1 <- length(lam1)
  n2 <- ifelse(is.null(dim(lam2)), length(lam2), dim(lam2)[2])
  n3 <- length(gamma)
  
  # The number of errors
  nerr <- 0
  code <- 0
  
  ## End of timing in the first stage.
  end_time <- Sys.time()
  time1 <- difftime(end_time, start_time, units = "secs")
  
  out_all <- lapply(1:nfolds, function(k){
    if(trace.it) cat(sprintf("Fold: %d/%d\n", k, nfolds))
    x_val <- x[fold==k,,drop=FALSE]
    y_val <- y[fold==k]
    
    x_train <- x[fold!=k,,drop=FALSE]
    y_train <- y[fold!=k]
    yclass_train <- yclass[fold!=k]
    
    if(is.null(M_fold) || is.null(U_fold)){
      # Fit with robustSIR function
      fit_fold <- robustSIR(x_train, y_train, yclass = yclass_train, lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, categorical = categorical, ...)
    }
    else{
      fit_fold <- robustSIR(M = M_fold[[k]], U = U_fold[[k]], nobs = sum(fold!=k), lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, categorical = categorical, ...)
    }
    
    err <- 0
    if(fit_fold$code != 0){
      err <- 1
      if(warning.it >= 2) warning(paste0("Fold ", k, ": No converged results returned from the robustSIR function."))
      eval_fold <- rep(NA_real_, length(fit_fold$beta))
      out <- list(eval_fold, err)
      return(out)
    }
    
    beta_l <- fit_fold$beta
    rank_l <- fit_fold$rank
    step_l <- fit_fold$step
    time_l <- fit_fold$time
    
    eval_fold <- eval_dc(beta_l, x_val, y_val, H = H, categorical = categorical, original_y = original_y)
    ind <- which(sapply(beta_l, is.null))
    rank_l[ind] <- -1
    eval_fold[ind] <- max(eval_fold, na.rm = TRUE)
    
    if(plot){
      dat <- data.frame(x = 1:length(eval_fold), y = eval_fold, rank = as.factor(rank_l))
      g <- ggplot(dat, aes(x = x, y = y, col = rank, shape = rank))+
        geom_point(size = 2)+
        scale_shape_manual(values = c(1:3,1:3,1:3))+
        labs(title=paste0("Fold ", k), x="", y="CV error")+
        theme_bw()+
        theme(legend.position = "bottom",
              legend.title = element_blank(),
              legend.text=element_text(size = 12),
              legend.key.width = unit(1, 'cm'),
              plot.title = element_text(size = 16, hjust = 0.5),
              axis.text = element_text(size = 16),
              axis.title = element_text(size = 16))
      print(g)
    }
    out <- list(eval_fold, err)
    out
  })
  
  eval_all <- do.call(rbind, lapply(out_all, "[[", 1))
  errs <- do.call(c, lapply(out_all, "[[", 2))
  nerr <- sum(errs)
  
  if((nerr != 0) && (nerr != nfolds)){
    code <- 3
    if(warning.it >= 2) warning(paste0("No converged results returned in", nerr, "folds."))
  }
  else if(nerr == nfolds){
    code <- 4
    if(warning.it >= 2) warning("No converged results returned in any fold.")
    return(list(beta = NULL, code = code))
  }
  
  if(is.vector(eval_all)){
    eval_all <- as.matrix(eval_all)
  }
  
  cvm <- colMeans(eval_all, na.rm=TRUE)
  cvsd <- sqrt(colMeans(scale(eval_all, cvm, FALSE)^2, na.rm = TRUE)/(nfolds-1))
  # min
  id_min <- which.min(cvm)
  id_lam1_min <- ceiling(id_min/(n2*n3))
  id_gamma_min <- ceiling((id_min-(id_lam1_min-1)*(n2*n3))/n2)
  id_lam2_min <- id_min-(id_lam1_min-1)*(n2*n3)-(id_gamma_min-1)*n2
  # lam1_min <- lam1[id_lam1_min]
  # gamma_min <- gamma[id_gamma_min]
  # lam2_min <- ifelse(is.null(dim(lam2)), lam2[id_lam2_min], lam2[id_gamma_min,id_lam2_min])
  
  ## Solution path
  if(solution.path || lambda.type == "1se.rank" || lambda.type == "1se.s"){
    if(is.null(M) || is.null(U)){
      # Fit with robustSIR function
      fit <- robustSIR(x, y, yclass = yclass, lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, categorical = categorical, true.beta = true.beta, ...)
    }
    else{
      fit <- robustSIR(M = M, U = U, nobs = NROW(x), lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, categorical = categorical, true.beta = true.beta, ...)
    }
    beta_path <- fit$beta
    rank_path <- fit$rank
    s_path <- fit$s
    dist_path <- fit$dist
    
    # lambda.type == "1se.rank"
    id_1se.rank <- which(cvm <= cvm[id_min] + cvsd[id_min])
    id_1se.rank <- id_1se.rank[!is.na(rank_path[id_1se.rank]) & !is.na(s_path[id_1se.rank])]
    id_1se.rank <- id_1se.rank[which(rank_path[id_1se.rank] == min(rank_path[id_1se.rank]))]
    id_1se.rank <- id_1se.rank[which(s_path[id_1se.rank] == min(s_path[id_1se.rank]))]
    id_1se.rank <- id_1se.rank[1]
    id_lam1_1se.rank <- ceiling(id_1se.rank/(n2*n3))
    id_gamma_1se.rank <- ceiling((id_1se.rank-(id_lam1_1se.rank-1)*(n2*n3))/n2)
    id_lam2_1se.rank <- id_1se.rank-(id_lam1_1se.rank-1)*(n2*n3)-(id_gamma_1se.rank-1)*n2
    
    # lambda.type == "1se.s"
    id_1se.s <- which(cvm <= cvm[id_min] + cvsd[id_min])
    id_1se.s <- id_1se.s[!is.na(rank_path[id_1se.s]) & !is.na(s_path[id_1se.s])]
    id_1se.s <- id_1se.s[which(s_path[id_1se.s] == min(s_path[id_1se.s]))]
    id_1se.s <- id_1se.s[which(rank_path[id_1se.s] == min(rank_path[id_1se.s]))]
    id_1se.s <- id_1se.s[1]
    id_lam1_1se.s <- ceiling(id_1se.s/(n2*n3))
    id_gamma_1se.s <- ceiling((id_1se.s-(id_lam1_1se.s-1)*(n2*n3))/n2)
    id_lam2_1se.s <- id_1se.s-(id_lam1_1se.s-1)*(n2*n3)-(id_gamma_1se.s-1)*n2
    
    ind <- which(sapply(rank_path, is.na))
    rank_path[ind] <- -1
    if(plot){
      dat <- data.frame(x = 1:length(cvm), y = cvm, rank = as.factor(rank_path))
      g <- ggplot(dat, aes(x = x, y = y, col = rank, shape = rank))+
          geom_point(size = 2)+
          scale_shape_manual(values = c(1:3,1:3,1:3))+
          labs(title="Solution path: overall CV error", x="", y="CV error")+
          theme_bw()+
          theme(legend.position = "bottom",
                legend.title = element_blank(),
                legend.text=element_text(size = 12),
                legend.key.width = unit(1, 'cm'),
                plot.title = element_text(size = 16, hjust = 0.5),
                axis.text = element_text(size = 16),
                axis.title = element_text(size = 16))
      print(g)
      if(solution.path.distance){
        dat <- data.frame(x = 1:length(cvm), y = dist_path, rank = as.factor(rank_path))
        g <- ggplot(dat, aes(x = x, y = y, col = rank, shape = rank))+
          geom_point(size = 2)+
          scale_shape_manual(values = c(1:3,1:3,1:3))+
          labs(title="Solution path: subspace distance", x="", y="Subspace distance")+
          theme_bw()+
          theme(legend.position = "bottom",
                legend.title = element_blank(),
                legend.text=element_text(size = 12),
                legend.key.width = unit(1, 'cm'),
                plot.title = element_text(size = 16, hjust = 0.5),
                axis.text = element_text(size = 16),
                axis.title = element_text(size = 16))
        print(g)
      }
      if(solution.path.s){
        dat <- data.frame(x = 1:length(cvm), y = s_path, rank = as.factor(rank_path))
        g <- ggplot(dat, aes(x = x, y = y, col = rank, shape = rank))+
          geom_point(size = 2)+
          scale_shape_manual(values = c(1:3,1:3,1:3))+
          labs(title="Solution path: sparsity level", x="", y="Sparsity level")+
          theme_bw()+
          theme(legend.position = "bottom",
                legend.title = element_blank(),
                legend.text=element_text(size = 12),
                legend.key.width = unit(1, 'cm'),
                plot.title = element_text(size = 16, hjust = 0.5),
                axis.text = element_text(size = 16),
                axis.title = element_text(size = 16))
        print(g)
      }
     }
  }else{
    beta_path <- NA
    rank_path <- NA
    s_path <- NA
    dist_path <- NA
  }
  
  if(lambda.type == "min"){
    id_lam1 <- id_lam1_min
    id_lam2 <- id_lam2_min
    id_gamma <- id_gamma_min
  }else if(lambda.type == "1se.rank"){
    id_lam1 <- id_lam1_1se.rank
    id_lam2 <- id_lam2_1se.rank
    id_gamma <- id_gamma_1se.rank
  }else if(lambda.type == "1se.s"){
    id_lam1 <- id_lam1_1se.s
    id_lam2 <- id_lam2_1se.s
    id_gamma <- id_gamma_1se.s
  }
  
  lam1_hat <- lam1[id_lam1]
  gamma_hat <- gamma[id_gamma]
  lam2_hat <- ifelse(is.null(dim(lam2)), lam2[id_lam2], lam2[id_gamma,id_lam2])
  
  start_time <- Sys.time() # Record time
  
  # Fit with robustSIR function
  if(is.null(M) || is.null(U)){
    fit <- robustSIR(x, y, yclass = yclass, lam1 = lam1_hat, lam2 = lam2_hat, gamma = gamma_hat, eps = eps, maxit = maxit, d = d, warning.it =  warning.it, categorical = categorical, ...)
  }
  else{
    fit <- robustSIR(M = M, U = U, nobs = NROW(x), lam1 = lam1_hat, lam2 = lam2_hat, gamma = gamma_hat, eps = eps, maxit = maxit, d = d, warning.it =  warning.it, categorical = categorical, ...)
  }
  
  if(fit$code != 0){
    code <- 5
    if(warning.it >= 1) warning("The estimated beta is null.")
    return(list(beta = NULL, code = code))
  }
  
  B <- fit$B
  beta <- fit$beta
  rank <- fit$rank
  
  ## End of timing in the second stage.
  end_time <- Sys.time()
  time2 <- difftime(end_time, start_time, units = "secs")

  ## The overall time.  
  time <- time1 + time2
  
  output <- list(beta = beta, B = B, rank = rank, eval = eval_all, cvm = cvm, id_lam1=id_lam1, id_lam2 = id_lam2, id_gamma = id_gamma, lam1 = lam1, lam2 = lam2, gamma = gamma, lam1_hat = lam1_hat, lam2_hat = lam2_hat, gamma_hat = gamma_hat, code = code, time = time, fold = fold, beta_path = beta_path, rank_path = rank_path, s_path = s_path, gammahat = gammahat, what = what)
  output
}


# ------------------ revised "msda" function from the "msda" package ---------------------- #
# * Following the idea of Zeng et al. (2024), which generates the sequence of tuning parameter candidates for their SEAS algorithm via the "msda" function, we also generates the candidate sequences for lambda1, lambda2 and the ADMM parameter gamma via the "msda" function.
# * We tailor the original "msda" function to adapt to our optimization problem in eq. (10) by replacing the kernel matrices "sigma" and "delta" in the original "msda" function with the Sigma and U matrices in eq. (9).
# * The inputs of the revised "msda" function are similar to those in its original version, thus we only describe the distinct ones. Please refer to the help document of "msda::msda" for details of other arguments.
#
# * Inputs:
# ========
# H: The number of slices.
# yclass, categorical: The arguments passed the "MU" function.
# M, U: The user-defined matrices, corresponding to the Sigma and U matrices in eq. (9).
# nobs: The number of observations.
# nclass: The number of classes in "y" if "categorical=TRUE" and equals to "H" otherwise.
# warning.it: The warning level. The higher this value is, the more warning messages will be printed.
# 
# * Outputs:
# ========
# outlist: The output list. Please refer to the "formatoutput" function in the "msda" package.
# x, y: The original predictor and response data.
# npasses: The number of total iterations.
# jerr: The error code from the "msda" Fortran code.
# M, U: The Sigma and U matrices returned from the "MU" function.
# rank: The list of ranks of each estimated matrix from MSDA.
# gammahat: Estimated gamma in eq (5).
# what: Estimate W in eq (7).

msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, lambda.factor=NULL, nlambda=NULL, lambda=NULL, dfmax=NULL, pmax=NULL, pf=NULL, M = NULL, U = NULL, nobs=NULL, nclass=NULL, eps=1e-04, maxit=1e+06, sml=1e-06, warning.it=0, perturb=NULL){
  ## Be careful about the value of lambda.factor since we add nobs = NULL in argument list.
  if(is.null(M) || is.null(U)){
    if(missing(x) || missing(y)) stop("Missing x or y.")
    if(is.data.frame(x)) x <- as.matrix(x)
    if(is.null(yclass)){
      if(categorical == FALSE){
        ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
        yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      }
      else if(categorical == TRUE){
        yclass <- y
      }
    }
    # if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    nclass <- as.integer(length(unique(yclass)))
    MU_out <- MU(x, y, yclass, categorical = categorical)
    M <- MU_out$M
    U <- MU_out$U
    what <- MU_out$what
    gammahat <- MU_out$gammahat
    nobs <- as.integer(dim(x)[1])
    nvars <- as.integer(dim(x)[2])
  }
  else{
    if(is.null(nobs)) stop("Missing nobs.")
    if(is.null(nclass)) stop("Missing nclass.")
    nvars <- NCOL(M)
  }
  
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.5, 1e-03)
  if(is.null(nlambda)) nlambda <- 100
  if(is.null(dfmax)) dfmax <- nobs
  if(is.null(pmax)) pmax <- min(dfmax*2 + 20, nvars)
  if(is.null(pf)) pf <- rep(1, nvars)
  if (!is.null(perturb)) 
    diag(M) <- diag(M) + perturb
  if(warning.it >= 4){
    verbose <- TRUE
  }else{
    verbose <- FALSE
  }
  H <- as.integer(dim(U)[2])
  ## parameter setup
  if (length(pf) != nvars) 
    stop("The size of penalty factor must be same as the number of input variables")
  maxit <- as.integer(maxit)
  verbose <- as.integer(verbose)
  sml <- as.double(sml)
  pf <- as.double(pf)
  eps <- as.double(eps)
  dfmax <- as.integer(dfmax)
  pmax <- as.integer(pmax)
  ## lambda setup
  nlam <- as.integer(nlambda)
  if (is.null(lambda)) {
    if (lambda.factor >= 1)
      stop("lambda.factor should be less than 1")
    flmin <- as.double(lambda.factor)
    ulam <- double(1)  #ulam=0 if lambda is missing
  } else {
    # flmin=1 if user define lambda
    flmin <- as.double(1)
    if (any(lambda < 0))
      stop("lambdas should be non-negative")
    ulam <- as.double(rev(sort(lambda)))  #lambda is declining
    nlam <- as.integer(length(lambda))
  }
  ## call Fortran core
  fit <- .Fortran("msda", obj = double(nlam), H, nvars, as.double(M), as.double(t(U)), pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam = integer(1), theta = double(pmax * H * nlam), itheta = integer(pmax), ntheta = integer(nlam),alam = double(nlam), npass = integer(1), jerr = integer(1))
  ## output
  outlist <- formatoutput(fit, maxit, pmax, nvars, H, warning.it)
  rank <- rep(NA_integer_, length(outlist$theta))
  ##
  for (i in 1:length(outlist$theta)){
    if(!is.null(outlist$theta[[i]])){
      rank[i] <- rank_func(outlist$theta[[i]], thrd = 1e-3)
    }
  }
  if(is.null(lambda))
    outlist$lambda <- lamfix(outlist$lambda)
  outlist <- c(outlist, list(x = x, y = y, npasses = fit$npass, jerr = fit$jerr, M = M, U = U, rank = rank, gammahat = gammahat, what = what))
  class(outlist) <- c("msda")
  outlist
}

# ------------------ Cross-validation function for the "msda" ---------------------- #
# * The inputs and outputs are similar to those in "msda" function. Only the distinct function are described in the following. Please refer to the description of the "msda" function for other input/output arguments.
#
# * Inputs:
# ========
# nfolds: The number of folds in the cross-validation.
# fold: User-defined fold index.
# warning.it: The warning level. The higher this value is, the more warning messages will be printed. This argument is passed to the "msda" function.
# plot: A logical argument. If TRUE, plots a cv-error curve.
# original_y: A logical value. If TRUE, use the original y in computing the distance correlation, which is used as the tuning parameter selection criterion. Otherwise, use the discretized y in computing the distance correlation. This argument is passed to the "eval_dc" function.
#
# * Outputs:
# ========
# beta: The optimal estimated matrix.
# id: The index of the optimal tuning parameter.
# lambda: The lambda sequence.
# lam_min: The optimal tuning parameter.
# rank: The rank of the optimal estimated matrix.
# M, U: The Sigma and U matrices obtained from the full data.
# M_fold, U_fold: Two lists containing the Sigma and U matrices obtained from each training fold.

cv.msda <- function(x, y, yclass=NULL, categorical=FALSE, H=5, lambda.factor=NULL, nlambda=100, nfolds=5, fold = NULL, lambda=NULL, maxit=1e3, warning.it = 1, plot=FALSE, original_y = FALSE){
  if(is.data.frame(x)) x <- as.matrix(x)
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  if(any(table(yclass) < nfolds)){
    stop(sprintf("The sample size of class %d is less than nfolds\n", which(table(yclass) < 5)))
  }
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  nclass <- length(unique(yclass))
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.5, 1e-03)
  # Fit full data, obtain the msda lambda candidates
  fit <- msda(x, y, yclass = yclass, lambda.factor = lambda.factor, nlambda = nlambda, lambda = lambda, maxit=maxit, warning.it = warning.it, categorical = categorical)
  gammahat <- fit$gammahat
  lambda <- fit$lambda
  beta_l <- fit$theta
  M <- fit$M
  U <- fit$U
  what <- fit$what
  rank_l <- fit$rank
  beta_l <- cut_mat(beta_l, 1e-3, rank_l)
  
  # Cross-validation
  if(is.null(fold)){
    ord <- order(y)
    y <- y[ord]
    yclass <- yclass[ord]
    x <- x[ord,]
    count <- as.numeric(table(yclass))
    fold <- c()
    for(cnt in count){
      fold <- c(fold, sample(rep(seq(nfolds), length = cnt)))
    }
  }
  else{
    nfolds <- length(unique(fold))
  }
  
  cv_out <- lapply(1:nfolds, function(k){
    x_train <- x[fold!=k,,drop=FALSE]		
    x_val <- x[fold==k,,drop=FALSE]		
    y_train <- y[fold!=k]
    y_val <- y[fold==k]
    yclass_train <- yclass[fold!=k]
    
    # matrix is already cut inside msda_func
    fit_fold <- msda(x_train, y_train, yclass_train, lambda.factor=lambda.factor, nlambda=nlambda, lambda = lambda, maxit=maxit, warning.it = warning.it, categorical = categorical)
    M_fold <- fit_fold$M
    U_fold <- fit_fold$U
    beta_fold <- fit_fold$theta
    rank_fold <- fit_fold$rank
    beta_fold <- cut_mat(beta_fold, 1e-3, rank_fold)
    
    # return evaluation of each fold
    eval_fold <- eval_dc(beta_fold, x_val, y_val, H = H, categorical = categorical, original_y = original_y)
    if(length(eval_fold) != length(lambda)){
      eval_fold <- c(eval_fold, rep(NA, length(lambda) - length(eval_fold)))
    }
    list(eval = eval_fold, M = M_fold, U = U_fold)
  })
  
  eval_all <- do.call(rbind, lapply(cv_out, "[[", 1))
  M_fold <- lapply(cv_out, "[[", 2)
  U_fold <- lapply(cv_out, "[[", 3)
  if(is.vector(eval_all)){
    eval_all <- t(as.matrix(eval_all))
  }
  
  ## No matrix is converged in any fold
  if(all(is.na(eval_all))) return(NULL)
  
  cvm <- colMeans(eval_all, na.rm = TRUE)
  # The optimal lambda1
  id_min <- which.min(cvm)
  lam_min <- lambda[id_min]
  beta <- as.matrix(beta_l[[id_min]])
  
  # Recalculate the rank
  rank <- rank_func(beta, thrd = 1e-3)
  
  if(plot){
    dat <- data.frame(x = 1:length(cvm), y = cvm)
    g <- ggplot(dat, aes(x = x, y = y))+
      geom_point(size = 2)+
      scale_shape_manual(values = c(1:3,1:3,1:3))+
      labs(x="", y="CV error")+
      theme_bw()+
      theme(legend.position = "bottom",
            legend.title = element_blank(),
            legend.text=element_text(size = 12),
            legend.key.width = unit(1, 'cm'),
            plot.title = element_text(size = 16, hjust = 0.5),
            axis.text = element_text(size = 16),
            axis.title = element_text(size = 16))
    print(g)
  }
  
  list(beta = beta, id = id_min, lambda = lambda, lam_min = lam_min, rank = rank, M = M, U = U, M_fold = M_fold, U_fold = U_fold, gammahat = gammahat, what = what)
}

# ---------- Compute the Sigma and U matrices in eq. (9) ------------- #
# * Input:
# x: n x p observation matrix for predictor.
# y: n-dimensional observation vector for response.
# yclass: Discretized response taking values in {1,...,H}.
# categorical: A logical value indicating whether y is categorical.
# H: The number of slices. The default value is 5.
#
# * Output:
# =========
# M, U: The estimated Sigma and U matrices in eq. (9).
# nclass: The number of classes in "y" if "categorical=TRUE" and equals to "H" otherwise.
# prior: The proportion of each class.
# gammahat: Estimated gamma in eq (5).
# what: Estimate W in eq (7).

MU <- function(x, y, yclass=NULL, categorical = FALSE, H = 5){
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
      nclass <- as.integer(length(unique(yclass)))
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  cls <- sort(unique(yclass))
  nclass <- length(cls)
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  prior <- sapply(cls, function(i){mean(yclass == i)})
  
  xbar <- colMeans(x)
  xc <- x - t(xbar)[rep(1, NROW(x)),]
  gamma <- apply(xc, 1, function(t){sqrt(mean(t^2))})  # gamma hat.
  weight <- 1/gamma^2
  mu <- c((t(weight) %*% x)/sum(weight))
  w <- (x - t(mu)[rep(1,NROW(x)),])/(gamma %*% t(rep(1,nvars)))
  
  wbar <- colMeans(w)
  wc <- (w - t(wbar)[rep(1,NROW(w)), ])
  
  M <- crossprod(wc)/nobs
  U <- matrix(0, nvars, nclass)
  for (i in 1:nclass){
    U[, i] <- colMeans(wc[yclass == cls[i],,drop = FALSE])
  }
  list(M = M, U = U, nclass = nclass, prior = prior, gammahat = gamma, what = w)
}

# --------------- ADMM algorithm for solving problem in eq. (10) ----------- #
# * This function is the main program for solving the optimization in eq. (10), the formulation of subspace estimation in ROSE method.
#
# * Inputs:
# =======
# M, U: The Sigma and U matrices in eq. (10).
# nobs: The number of observations.
# nvars: The number of predictors.
# lam1, lam2, gamma: The user-specified sequences of tuning parameter lambda1, lambda2, and the ADMM parameter gamma.
# eps: The tolerance of convergence in ADMM algorithm.
# maxit: The maximal iterations in ADMM algorithm.
# d: The true structural dimension. This argument is used when the oracle knowledge of the structural dimension is given.
# warning.it: The warning level. The higher this value is, the more warning messages will be printed. This argument is passed to the "msda" function.
# ...: Other arguments.
# 
# * Outputs:
# ========
# beta: A list of the estimated basis matrices of central subspace.
# B: A list of the estimated B matrices.
# rank: A vector of the estimated ranks.
# s: A vector of the estimated sparsity levels.
# step: A vector of the number of iterations corresponding to each tuning parameter set.
# time: The execution time of implementing the ADMM algorithm for solving eq (10) once.
# nlam: The number of converged matrices.
# sv_list_B, sv_list_C: Two lists of the singular values of the converged matrices B and C.

admm <- function(M, U, nobs, nvars, lam1, lam2, gam, eps=1e-3, maxit=1e+3, d = NULL, warning.it = 0, ...){
  # since the user is required to provide lam1, then set flmin=1
  opts <- list(...)
  if(is.null(opts$nlam)) opts$nlam <- as.integer(1)
  if(is.null(opts$H)) opts$H <- as.integer(dim(U)[2])
  if(is.null(opts$nvars)) opts$nvars <- as.integer(nvars)
  if(is.null(opts$pf)) opts$pf <- as.double(rep(1, nvars))
  if(is.null(opts$dfmax)) opts$dfmax <- as.integer(nobs)
  if(is.null(opts$pmax)) opts$pmax <- as.integer(min(nobs*2+20, nvars))
  if(is.null(opts$flmin)) opts$flmin <- as.double(1)
  if(is.null(opts$eps_inner)) opts$eps_inner <- as.double(1e-04)
  if(is.null(opts$maxit_inner)) opts$maxit_inner <- as.integer(1e+6)
  if(is.null(opts$sml)) opts$sml <- as.double(1e-6)
  if(is.null(opts$verbose)) opts$verbose <- as.integer(FALSE)
  if(is.null(opts$nalam)) opts$nalam <- integer(1)
  if(is.null(opts$theta)) opts$theta <- double(opts$pmax * opts$H * opts$nlam)
  if(is.null(opts$itheta)) opts$itheta <- integer(opts$pmax)
  if(is.null(opts$ntheta)) opts$ntheta <- integer(opts$nlam)
  if(is.null(opts$alam)) opts$alam <- double(opts$nlam)
  if(is.null(opts$npass)) opts$npass <- integer(1)
  if(is.null(opts$jerr)) opts$jerr <- integer(1)
  
  M0 <- M
  U0 <- U
  n1 <- length(lam1)
  n2 <- ifelse(is.null(dim(lam2)), length(lam2), ncol(lam2))
  n3 <- length(gam)
  nparams <- n1*n2*n3
  
  B_l <- vector("list", nparams)
  beta_l <- vector("list", nparams)
  step_l <- rep(NA_integer_, nparams) # To store the iteration times of each run
  time_l <- rep(NA_real_, nparams) # To store the running time of each run
  rank_l <- rep(NA_integer_, nparams)
  s_l <- rep(NA_integer_, nparams)
  
  sv_list_B <- vector("list", nparams)
  sv_list_C <- vector("list", nparams)
  
  # The number of converged matrices
  nlam_cvg <- 0
  
  for(i in 1:n1){
    lambda1 <- as.double(lam1[i])
    
    for(j in 1:n3){
      gamma <- gam[j]
      
      for(k in 1:n2){
        lambda2 <- ifelse(is.null(dim(lam2)), lam2[k], lam2[j,k])
        
        M <- M0 + gamma*diag(rep(1,ncol(M0)), ncol(M0),ncol(M0))
        
        # Initialize three matrices
        Bold <- matrix(0,dim(U0)[1], dim(U0)[2])
        Cold <- matrix(0,dim(U0)[1], dim(U0)[2])
        etaold <- matrix(0,dim(U0)[1], dim(U0)[2])
        
        # The MAIN loop of admm method
        step <- 0    
        start_time <- Sys.time()
        
        repeat{
          step <- step + 1
          
          # Update B
          U <- U0 - etaold + gamma * Cold
          out_B <- updateB(M, U, lambda1, opts)
          Bnew <- out_B$Bnew
          jerr <- out_B$jerr
          if(jerr != 0) break
          
          # Update C
          Cnew <- updateC(Bnew, lambda2, gamma, etaold)
          
          # Update U
          etanew <- etaold + gamma * (Bnew - Cnew)
          
          # Code 1: success
          if(max(abs(Bnew - Bold)) < eps){
            jerr <- 1
            break
          }
          # Code 404: then maximal iteration is reached
          if(step > maxit){
            jerr <- 404
            if(warning.it >= 3) warning('Maximal iteration is reached.')
            break
          }
          Bold <- Bnew
          Cold <- Cnew
          etaold <- etanew
        }# End of repeat 
        end_time <- Sys.time()  # The time for each repeat
        time <- difftime(end_time, start_time, units = "secs")
        # Code < -10000: non-sparse matrix
        if(jerr < -10000){
          if(warning.it >= 3) warning(paste0('The ',i, 'th lam1 is too small, no execution for the rest parameters lam2 and gamma.'))
          break
        }
        # Code 1: success, save the matrix and the related information.
        if(jerr==1){
          index <- (i-1)*n2*n3 + (j-1)*n2 + k
          nlam_cvg <- nlam_cvg + 1
          B_l[[index]] <- Bnew
          step_l[index] <- step
          time_l[index] <- time
          if(is.null(d)) rank <- rank_func(Cnew, thrd = 1e-3)
          else rank <- d
          rank_l[index] <- rank
          # save the singular values of each candidates matrix B and C
          sv_list_B[[index]] <- svd(Bnew)$d
          sv_list_C[[index]] <- svd(Cnew)$d
          # Cut and select the left singular vector of Bnew
          if(rank == 0){
            beta <- matrix(0, nrow(Bnew), ncol(Bnew))
          }else{
            tmp <- svd(Bnew)$u[,1:rank, drop = FALSE]
            vec <- as.vector(tmp)
            vec[abs(vec) < 1e-3] <- 0
            beta <- matrix(vec, nrow(tmp), ncol(tmp))
          }
          beta_l[[index]] <- beta
          var_ind <- apply(beta, 1, function(x){any(x!=0)})
          s_l[index] <- sum(var_ind)
        }
      }# End of lambda2
      if(jerr < -10000) break
    }# End of gam
  }# End of lambda1
  return(list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, step = step_l, time = time_l, nlam = nlam_cvg, sv_list_B = sv_list_B, sv_list_C = sv_list_C))
}

# ------------ Update B iterate in the "admm" algorithm ----------- #
# * Inputs:
# =======
# M, U: The Sigma and U matrices in eq. (10).
# lambda1: The user-specified sequences of tuning parameter lambda1.
# opts: The options passed from the "admm" algorithm.
#
# * Output:
# =======
# Bnew: The converged matrix.
# jerr: The error code from the "msda" Fortran code.

updateB <- function(M, U, lambda1, opts){
  U <- t(U)
  fit <- .Fortran("msda", obj = opts$nlam, opts$H, opts$nvars, as.double(M), as.double(U), opts$pf, opts$dfmax, opts$pmax, opts$nlam, opts$flmin, lambda1, opts$eps_inner, opts$maxit_inner, opts$sml, opts$verbose, nalam = opts$nalam, theta = opts$theta, itheta = opts$itheta, ntheta = opts$ntheta, alam = opts$alam, npass = opts$npass, jerr = opts$jerr)
  if(fit$jerr != 0){return(list(Bnew = NULL, jerr = fit$jerr))} # Code: non-zero, abnormal result.
  outlist <- formatoutput(fit, opts$maxit_inner, opts$pmax, opts$nvars, opts$H)
  Bnew <- as.matrix(outlist$theta[[1]])
  list(Bnew = Bnew, jerr = fit$jerr)
}

# ------------ Update C iterate in the "admm" function ----------- #
# * Inputs:
# =======
# Bnew: The B iterate updated from the "updateB" function.
# lambda2: The lambda2 tuning parameter.
# gamma: The ADMM tuning parameter.
# etaold: The dual parameter.
#
# * Output:
# =======
# Cnew: The converged matrix.

updateC <- function(Bnew, lambda2, gamma, etaold){
  Btemp <- Bnew + 1/gamma * etaold
  svd_B <- svd(Btemp)
  lamtemp <- pmax(0, svd_B$d-lambda2/gamma)
  Cnew <- svd_B$u %*% diag(lamtemp, nrow = length(lamtemp), ncol = length(lamtemp)) %*% t(svd_B$v)
  Cnew
}

