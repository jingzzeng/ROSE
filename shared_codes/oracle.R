oracle <- function(x = NULL, y = NULL, yclass = NULL, d = NULL, categorical=FALSE, H=5, M = NULL, U = NULL, nobs = NULL, nlambda=NULL, lambda.factor=NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.3, length.out = 10), lam2_fac=seq(0.01,0.3, length.out = 10), eps = 1e-3, maxit = 1e+3, warning.it = 1, ...){
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
    # if(any(table(yclass) < 5)){
    #   if(warning.it >= 1) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    # }
    if(is.null(gamma)){
      gamma <- c(10,30,50)
    }
    if(is.null(lam1) || is.null(lam2)){
      fit_1 <- cv.msda_oracle(x, y, yclass = yclass, nlambda=nlambda, lambda.factor=lambda.factor, nfolds = 5, maxit=1e3, warning.it = warning.it, categorical = categorical)
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
      MU_out <- MU_oracle(x, y, yclass)
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
    # fit <- admm_oracle(M, U, nobs, nvars, lam1, lam2, gamma, eps, maxit, d, warning.it = warning.it, ...)
    fit <- admm(M, U, nobs, nvars, lam1, lam2, gamma, eps, maxit, d, warning.it = warning.it, ...)
    B_l <- fit$B
    beta_l <- fit$beta
    if (all(sapply(beta_l, is.null))){
      code <- 1
      rank_l <- rep(NA, length(beta_l))
      s_l <- rep(NA, length(beta_l))
      if(warning.it >= 2) warning("In \"oracle\": No converged results returned.\n")
      if(length(B_l) == 1){
        B_l <- B_l[[1]]
        beta_l <- beta_l[[1]]
        rank_l <- NA
        s_l <- NA
      }
      return(list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, step = NA, time = NA, code = code))
    }
    rank_l <- fit$rank
    s_l <- fit$s
    step_l <- fit$step
    time_l <- fit$time
    if(length(B_l) == 1){
      B_l = B_l[[1]]; beta_l = beta_l[[1]]; rank_l = rank_l[[1]]; s_l = s_l[[1]]; step_l = step_l[[1]]; time_l = time_l[[1]]
    }
    output <- list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, lam1 = lam1, lam2 = lam2, gamma = gamma, step = step_l, time = time_l, code = code)
  }
  output
}


oracle_val <- function(x_train, y_train, x_val, y_val, yclass = NULL, d = NULL, categorical=FALSE, H=5, lambda.factor=NULL, nlambda=NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.3, length.out = 10), lam2_fac=seq(0.01,0.3, length.out = 10), plot = FALSE, eps = 1e-3, maxit = 1e+3, trace.it = TRUE, true_beta = NULL, warning.it = 1, ...){
  if(is.data.frame(x_train)) x_train <- as.numeric(x_train)
  if(is.data.frame(x_val)) x_val <- as.numeric(x_val)
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y_train, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y_train, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y_train
    }
  }
  # if(any(table(yclass) < 5)){
  #   if(warning.it >= 1) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
  # }
  
  ord <- order(y_train)
  y_train <- y_train[ord]
  yclass <- yclass[ord]
  x_train <- x_train[ord,]
  M <- U <- NULL
  
  if(is.null(gamma)){
    gamma <- c(10,30,50)
  }
  if(is.null(lam1) || is.null(lam2)){
    fit_1 <- msda_oracle(x_train, y_train, yclass = yclass, nlambda = nlambda, lambda.factor=lambda.factor, maxit=1e3, warning.it = warning.it)
    M <- fit_1$M
    U <- fit_1$U
    lam_msda <- fit_1$lambda
    beta_msda <- fit_1$theta
    rank_msda <- fit_1$rank
    # Cut matrix
    beta_msda <- cut_mat(beta_msda, 1e-3, rank_msda)
    eval_msda <- eval_dc_oracle(beta_msda, x_val, y_val, H = H, categorical = categorical)
    id_min_msda <- which.min(eval_msda)
    lam1_min_msda <- lam_msda[id_min_msda]
    beta_msda <- as.matrix(beta_msda[[id_min_msda]])
    if(plot){
      dat <- data.frame(x = 1:length(eval_msda), y = eval_msda)
      g <- ggplot(dat, aes(x = x, y = y))+
        geom_point(size = 2)+
        xlab("")+
        ylab("Distance correlation")+
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
    # lambda1 and lambda candidates
    if(is.null(lam1)) lam1 <- (lam1_min_msda)*lam1_fac
    if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
    # Code 1: msda matrix is zero matrix
    if (all(lam2 == 0)){
      lam2 <- 0
      if(warning.it >= 1) warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.")
    }
  }
  n1 <- length(lam1)
  n2 <- ifelse(is.null(dim(lam2)), length(lam2), dim(lam2)[2])
  n3 <- length(gamma)
  
  # Error code
  code <- 0
  
  # Fit with oracle function
  if (trace.it) cat("Training Fit\n")
  if(is.null(M) || is.null(U)){
    fit <- oracle(x_train, y_train, yclass, lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, ...)
  }
  else{
    fit <- oracle(M = M, U = U, nobs = NROW(x_train), lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, ...)
  }
  
  if(fit$code != 0){
    if(warning.it >= 1) warning("No converged results returned from the oracle function.")
    return(list(beta = NULL, code = fit$code))
  }
  
  B_l <- fit$B
  beta_l <- fit$beta
  rank_l <- fit$rank
  step_l <- fit$step
  time_l <- fit$time
  
  ## Solution Path
  dist_l <- rep(NA, length(beta_l))
  for(i in seq_along(beta_l)){
    if((rank_l[i] == 0) || (is.na(rank_l[i]))){dist_l[i] <- Inf}
    else{dist_l[i] <- subspace(true_beta, beta_l[[i]])}
  }
  id_path <- which.min(dist_l)
  beta_path <- beta_l[[id_path]]
  rank_path <- rank_l[id_path]
  dist_path <- dist_l[id_path]
  ##
  
  if (trace.it) cat("Validation\n")
  eval <- eval_dc_oracle(beta_l, x_val, y_val, H = H, categorical = categorical)
  ind <- which(sapply(beta_l, is.null))
  rank_l[ind] <- -1
  eval[ind] <- max(eval, na.rm = TRUE)
  
  if(plot){
    dat <- data.frame(x = 1:length(eval), y = eval, rank = as.factor(rank_l))
    g <- ggplot(dat, aes(x = x, y = y, col = rank,  shape = rank))+
      geom_point(size = 2)+
      xlab("")+
      ylab("Distance correlation")+
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
  
  # The optimal beta and tuning parameters
  id_min <- which.min(eval)
  B <- B_l[[id_min]]
  beta <- beta_l[[id_min]]
  rank <- rank_l[id_min]
  id_lam1 <- ceiling(id_min/(n2*n3))
  id_gamma <- ceiling((id_min-(id_lam1-1)*(n2*n3))/n2)
  id_lam2 <- id_min-(id_lam1-1)*(n2*n3)-(id_gamma-1)*n2
  lam1_min <- lam1[id_lam1]
  gamma_min <- gamma[id_gamma]
  lam2_min <- ifelse(is.null(dim(lam2)), lam2[id_lam2], lam2[id_gamma,id_lam2])
  
  ## Code 3: the estimated beta is null
  if(is.null(beta)){
    code <- 2
    if(warning.it >= 1) warning("The estimated beta is null.")
    return(list(beta = beta, code = code))
  }
  
  output <- list(beta = beta, B = B, rank = rank, eval = eval, M = M, U = U, id_lam1=id_lam1, id_lam2 = id_lam2, id_gamma = id_gamma, lam1 = lam1, lam2 = lam2, gamma = gamma, lam1_min = lam1_min, lam2_min = lam2_min, gamma_min = gamma_min, time = time_l, step = step_l, rank_l = rank_l, id_min = id_min, code = code, beta_path = beta_path, rank_path = rank_path)
  
  output
}

cv.oracle <- function(x, y, yclass = NULL, d = NULL, categorical=FALSE, H=5, lambda.factor=NULL, nlambda=NULL, nfolds = 5, fold = NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.3, length.out = 10), lam2_fac=seq(0.01,0.3, length.out = 10), plot = FALSE, eps = 1e-3, maxit = 1e+3, trace.it = TRUE, solution.path = FALSE, warning.it = 1, ...){
  
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
    fit_1 <- cv.msda_oracle(x, y, yclass = yclass, nlambda=nlambda, lambda.factor=lambda.factor, fold = fold, maxit=1e3, warning.it =  warning.it, plot = plot, categorical = categorical)
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
      # Fit with oracle function
      fit_fold <- oracle(x_train, y_train, yclass = yclass_train, lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, ...)
    }
    else{
      fit_fold <- oracle(M = M_fold[[k]], U = U_fold[[k]], nobs = sum(fold!=k), lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, ...)
    }
    
    err <- 0
    if(fit_fold$code != 0){
      err <- 1
      if(warning.it >= 2) warning(paste0("Fold ", k, ": No converged results returned from the oracle function."))
      eval_fold <- rep(NA_real_, length(fit_fold$beta))
      out <- list(eval_fold, err)
      return(out)
    }
    
    beta_l <- fit_fold$beta
    rank_l <- fit_fold$rank
    step_l <- fit_fold$step
    time_l <- fit_fold$time
    
    eval_fold <- eval_dc_oracle(beta_l, x_val, y_val, H = H, categorical = categorical)
    ind <- which(sapply(beta_l, is.null))
    rank_l[ind] <- -1
    eval_fold[ind] <- max(eval_fold, na.rm = TRUE)
    
    if(plot){
      dat <- data.frame(x = 1:length(eval_fold), y = eval_fold, rank = as.factor(rank_l))
      g <- ggplot(dat, aes(x = x, y = y, col = rank, shape = rank))+
        geom_point(size = 2)+
        labs(title=paste0("Fold ", k), x="", y="Distance correlation")+
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
  id_min <- which.min(cvm)
  id_lam1 <- ceiling(id_min/(n2*n3))
  id_gamma <- ceiling((id_min-(id_lam1-1)*(n2*n3))/n2)
  id_lam2 <- id_min-(id_lam1-1)*(n2*n3)-(id_gamma-1)*n2
  lam1_min <- lam1[id_lam1]
  gamma_min <- gamma[id_gamma]
  lam2_min <- ifelse(is.null(dim(lam2)), lam2[id_lam2], lam2[id_gamma,id_lam2])
  
  ## Solution path
  if(solution.path){
    if(is.null(M) || is.null(U)){
      # Fit with oracle function
      fit <- oracle(x, y, yclass = yclass, lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, categorical = categorical, ...)
    }
    else{
      fit <- oracle(M = M, U = U, nobs = NROW(x), lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, categorical = categorical, ...)
    }
    # beta_l <- fit$beta
    # rank_l <- fit$rank
    # s_l <- fit$s
    beta_path <- fit$beta
    rank_path <- fit$rank
    s_path <- fit$s
    
    # dist_l <- rep(NA, length(beta_l))
    # for(i in seq_along(beta_l)){
    #   if((rank_l[i] == 0) || (is.na(rank_l[i]))){dist_l[i] <- Inf}
    #   else{dist_l[i] <- subspace(true_beta, beta_l[[i]])}
    # }
    # id_path <- which.min(dist_l)
    # beta_path <- beta_l[[id_path]]
    # rank_path <- rank_l[id_path]
    # s_path <- s_l[id_path]
    # dist_path <- dist_l[id_path]
  }else{
    beta_path <- NA
    rank_path <- NA
    s_path <- NA
    # dist_path <- NA
  }
  
  # ind <- which(sapply(beta_l, is.null))
  # rank_l[ind] <- -1
  # if(plot){
  #   dat <- data.frame(x = 1:length(cvm), y = cvm, rank = as.factor(rank_l))
  #   g <- ggplot(dat, aes(x = x, y = y, col = rank, shape = rank))+
  #     geom_point(size = 2)+
  #     labs(title="Overall", x="", y="Distance correlation")+
  #     theme_bw()+
  #     theme(legend.position = "bottom",
  #           legend.title = element_blank(),
  #           legend.text=element_text(size = 12),
  #           legend.key.width = unit(1, 'cm'),
  #           plot.title = element_text(size = 16, hjust = 0.5),
  #           axis.text = element_text(size = 16),
  #           axis.title = element_text(size = 16))
  #   print(g)
  # }
  
  start_time <- Sys.time() # Record the time
  
  # Fit with oracle function
  if(is.null(M) || is.null(U)){
    fit <- oracle(x, y, yclass = yclass, lam1 = lam1_min, lam2 = lam2_min, gamma = gamma_min, eps = eps, maxit = maxit, d = d, warning.it =  warning.it, categorical = categorical, ...)
  }
  else{
    fit <- oracle(M = M, U = U, nobs = NROW(x), lam1 = lam1_min, lam2 = lam2_min, gamma = gamma_min, eps = eps, maxit = maxit, d = d, warning.it =  warning.it, categorical = categorical, ...)
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
  
  output <- list(beta = beta, B = B, rank = rank, eval = eval_all, id_lam1=id_lam1, id_lam2 = id_lam2, id_gamma = id_gamma, lam1 = lam1, lam2 = lam2, gamma = gamma, lam1_min = lam1_min, lam2_min = lam2_min, gamma_min = gamma_min, code = code, time = time, fold = fold, beta_path = beta_path, rank_path = rank_path, s_path = s_path)
  output
}

oracle_path <- function(x, y, yclass = NULL, d = NULL, categorical=FALSE, H=5, lambda.factor=NULL, nlambda=NULL, lam1 = NULL, lam2 = NULL, gamma = NULL, lam1_fac=seq(1.0,0.3, length.out = 10), lam2_fac=seq(0.01,0.3, length.out = 10), plot = FALSE, eps = 1e-3, maxit = 1e+3, trace.it = TRUE, true_beta = NULL, ...){
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
  
  ord <- order(y)
  y <- y[ord]
  yclass <- yclass[ord]
  x <- x[ord,]
  M <- U <- NULL
  
  nobs <- dim(x)[1]
  count <- as.numeric(table(yclass))
  
  if(is.null(gamma)){
    gamma <- c(10,30,50)
  }
  if(is.null(lam1) || is.null(lam2)){
    fit_1 <- msda_oracle_path(x, y, yclass = yclass, nlambda=nlambda, lambda.factor=lambda.factor, maxit=1e3, true_beta = true_beta, plot = plot)
    M <- fit_1$M
    U <- fit_1$U
    id_min_msda <- fit_1$id_path
    lam1_min_msda <- fit_1$lam_path
    beta_msda <- as.matrix(fit_1$beta_path)
    if(is.null(lam1)) lam1 <- (lam1_min_msda)*lam1_fac
    if(is.null(lam2)) lam2 <- svd(beta_msda)$d[1] * matrix(gamma, ncol = 1) %*% matrix(lam2_fac, nrow = 1)
    if (all(lam2 == 0)){
      warning("The automatically generated lambda 2 is zero, no nuclear norm penalty is imposed.")
    }
  }
  n1 <- length(lam1)
  n2 <- ifelse(is.null(dim(lam2)), length(lam2), dim(lam2)[2])
  n3 <- length(gamma)
  
  # The number of errors
  nerr <- 0
  code <- 0
  
  ##
  if(is.null(M) || is.null(U)){
    # Fit with oracle function
    fit <- oracle(x, y, yclass = yclass, lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, ...)
  }
  else{
    fit <- oracle(M = M, U = U, nobs = NROW(x), lam1 = lam1, lam2 = lam2, gamma = gamma, eps = eps, maxit = maxit, d = d, ...)
  }
  if(fit$code == 1){return(list(beta_path = NULL))}
  
  if(!is.list(fit$beta)) beta_l <- list(fit$beta)
  else beta_l <- fit$beta
  
  rank_l <- fit$rank
  s_l <- fit$s
  dist_l <- rep(NA, length(beta_l))
  for(i in seq_along(beta_l)){
    if((rank_l[i] == 0) || (is.na(rank_l[i]))){dist_l[i] <- Inf}
    else{dist_l[i] <- subspace(true_beta, beta_l[[i]])}
  }
  
  if(plot){
    dat <- data.frame(x = 1:length(dist_l), y = dist_l, rank = as.factor(rank_l))
    g <- ggplot(dat, aes(x = x, y = y, col = rank))+
      geom_point(size = 1)+
      labs(title="oracle", x="", y="Subspace distance")+
      theme_bw()+
      theme(axis.title = element_text(size=16),
            axis.text = element_text(size=14),
            plot.title = element_text(size=16))
    print(g)
  }
  
  id_path <- which.min(dist_l)
  id_lam1 <- ceiling(id_path/(n2*n3))
  id_gamma <- ceiling((id_path - (id_lam1-1)*(n2*n3))/n2)
  id_lam2 <- id_path - (id_lam1-1)*(n2*n3) - (id_gamma-1)*n2
  lam1_path <- lam1[id_lam1]
  lam2_path <- lam2[id_lam2]
  gamma_path <- gamma[id_gamma]
  beta_path <- beta_l[[id_path]]
  rank_path <- rank_l[id_path]
  
  output <- list(beta_path = beta_path, rank_path = rank_path, lam1 = lam1, lam2 = lam2, gamma = gamma, lam1_path = lam1_path, lam2_path = lam2_path, gamma_path = gamma_path, id_lam1 = id_lam1, id_lam2 = id_lam2, id_gamma = id_gamma)
  output
}


# calculate the Beta and the lambda sequences using msda function
msda_oracle <- function(x, y, yclass=NULL, categorical=FALSE, H=5, lambda.factor=NULL, nlambda=NULL, lambda=NULL, dfmax=NULL, pmax=NULL, pf=NULL, M = NULL, U = NULL, nobs=NULL, nclass=NULL, eps=1e-04, maxit=1e+06, sml=1e-06, warning.it=0, perturb=NULL){
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
    if(any(table(yclass) < 5)) warning(sprintf("The sample size of class %d is less than 5\n", which(table(yclass) < 5)))
    nclass <- as.integer(length(unique(yclass)))
    MU_out <- MU_oracle(x, y, yclass)
    M <- MU_out$M
    U <- MU_out$U
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
  outlist <- c(outlist, list(x = x, y = y, npasses = fit$npass, jerr = fit$jerr, M = M, U = U, rank = rank))
  class(outlist) <- c("msda")
  outlist
}

cv.msda_oracle <- function(x, y, yclass=NULL, categorical=FALSE, H=5, lambda.factor=NULL, nlambda=100, nfolds=5, fold = NULL, lambda=NULL, maxit=1e3, warning.it = 1, plot=FALSE){
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
  fit <- msda_oracle(x, y, yclass = yclass, lambda.factor = lambda.factor, nlambda = nlambda, lambda = lambda, maxit=maxit, warning.it = warning.it, categorical = categorical)
  lambda <- fit$lambda
  beta_l <- fit$theta
  M <- fit$M
  U <- fit$U
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
    fit_fold <- msda_oracle(x_train, y_train, yclass_train, lambda.factor=lambda.factor, nlambda=nlambda, lambda = lambda, maxit=maxit, warning.it = warning.it, categorical = categorical)
    M_fold <- fit_fold$M
    U_fold <- fit_fold$U
    beta_fold <- fit_fold$theta
    rank_fold <- fit_fold$rank
    beta_fold <- cut_mat(beta_fold, 1e-3, rank_fold)
    
    # return evaluation of each fold
    eval_fold <- eval_dc_oracle(beta_fold, x_val, y_val, H = H, categorical = categorical)
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
      xlab("")+
      ylab("Distance correlation")+
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
  
  list(beta = beta, id = id_min, lambda = lambda, lam_min = lam_min, rank = rank, M = M, U = U, M_fold = M_fold, U_fold = U_fold)
}


msda_oracle_path <- function(x, y, yclass=NULL, categorical=FALSE, H=5, lambda.factor=NULL, nlambda=100, lambda=NULL, maxit=1e3, true_beta = NULL, plot = FALSE){
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
  nobs <- as.integer(dim(x)[1])
  nvars <- as.integer(dim(x)[2])
  nclass <- length(unique(yclass))
  if(is.null(lambda.factor)) lambda.factor <- ifelse((nobs - nclass)<=nvars, 0.5, 1e-03)
  # Fit full data, obtain the msda lambda candidates
  fit <- msda_oracle(x, y, yclass = yclass, lambda.factor = lambda.factor, nlambda = nlambda, lambda = lambda, maxit=maxit)
  lambda <- fit$lambda
  beta_l <- fit$theta
  rank_l <- fit$rank
  beta_l <- cut_mat(beta_l, 1e-3, rank_l)
  M <- fit$M
  U <- fit$U
  
  # The best solution on the solution path.
  dist_l <- rep(NA, length(beta_l))
  for(i in seq_along(beta_l)){
    if((rank_l[i] == 0) || (is.na(rank_l[i]))){dist_l[i] <- Inf}
    else{dist_l[i] <- subspace(true_beta, beta_l[[i]])}
  }
  
  if(plot){
    dat <- data.frame(x = 1:length(dist_l), y = dist_l)
    g <- ggplot(dat, aes(x = x, y = y))+
      geom_point(size = 1)+
      labs(title="MSDA", x="", y="Subspace distance")+
      theme_bw()+
      theme(axis.title = element_text(size=16),
            axis.text = element_text(size=14),
            plot.title = element_text(size=16))
    print(g)
  }
  
  id_path <- which.min(dist_l)
  beta_path <- beta_l[[id_path]]
  rank_path <- rank_l[id_path]
  lam_path <- lambda[id_path]
  
  list(beta_path = beta_path, id_path = id_path, lambda = lambda, lam_path = lam_path, rank_path = rank_path, M = M, U = U)
}

MU_oracle <- function(x, y, yclass=NULL, categorical = FALSE, H = 5){
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
  # x_c <- x - matrix(xbar, nobs, nvars, byrow = TRUE)
  xc <- (x - t(xbar)[rep(1,NROW(x)), ])
  # M <- crossprod(x_c/sqrt(nobs))
  M <- crossprod(xc)/nobs
  
  U <- matrix(0, nvars, nclass)
  for (i in 1:nclass){
    U[, i] <- colMeans(xc[yclass == cls[i],, drop=FALSE])
  }
  list(M = M, U = U, nclass = nclass, prior=prior)
}


# # admm algorithm function
# admm_oracle <- function(M, U, nobs, nvars, lam1, lam2, gam, eps=1e-3, maxit=1e+3, d = NULL, warning.it = 0, ...){
#   # since the user is required to provide lam1, then set flmin=1
#   opts <- list(...)
#   if(is.null(opts$nlam)) opts$nlam <- as.integer(1)
#   if(is.null(opts$H)) opts$H <- as.integer(dim(U)[2])
#   if(is.null(opts$nvars)) opts$nvars <- as.integer(nvars)
#   if(is.null(opts$pf)) opts$pf <- as.double(rep(1, nvars))
#   if(is.null(opts$dfmax)) opts$dfmax <- as.integer(nobs)
#   if(is.null(opts$pmax)) opts$pmax <- as.integer(min(nobs*2+20, nvars))
#   if(is.null(opts$flmin)) opts$flmin <- as.double(1)
#   if(is.null(opts$eps_inner)) opts$eps_inner <- as.double(1e-04)
#   if(is.null(opts$maxit_inner)) opts$maxit_inner <- as.integer(1e+6)
#   if(is.null(opts$sml)) opts$sml <- as.double(1e-6)
#   if(is.null(opts$verbose)) opts$verbose <- as.integer(FALSE)
#   if(is.null(opts$nalam)) opts$nalam <- integer(1)
#   if(is.null(opts$theta)) opts$theta <- double(opts$pmax * opts$H * opts$nlam)
#   if(is.null(opts$itheta)) opts$itheta <- integer(opts$pmax)
#   if(is.null(opts$ntheta)) opts$ntheta <- integer(opts$nlam)
#   if(is.null(opts$alam)) opts$alam <- double(opts$nlam)
#   if(is.null(opts$npass)) opts$npass <- integer(1)
#   if(is.null(opts$jerr)) opts$jerr <- integer(1)
#   
#   M0 <- M
#   U0 <- U
#   n1 <- length(lam1)
#   n2 <- ifelse(is.null(dim(lam2)), length(lam2), ncol(lam2))
#   n3 <- length(gam)
#   nparams <- n1*n2*n3
#   
#   B_l <- vector("list", nparams)
#   beta_l <- vector("list", nparams)
#   step_l <- rep(NA_integer_, nparams) # To store the iteration times of each run
#   time_l <- rep(NA_real_, nparams) # To store the running time of each run
#   rank_l <- rep(NA_integer_, nparams)
#   s_l <- rep(NA_integer_, nparams)
#   
#   sv_list_B <- vector("list", nparams)
#   sv_list_C <- vector("list", nparams)
#   
#   # The number of converged matrices
#   nlam_cvg <- 0
#   
#   for(i in 1:n1){
#     lambda1 <- as.double(lam1[i])
#     
#     for(j in 1:n3){
#       gamma <- gam[j]
#       
#       for(k in 1:n2){
#         lambda2 <- ifelse(is.null(dim(lam2)), lam2[k], lam2[j,k])
#         
#         M <- M0 + gamma*diag(rep(1,ncol(M0)), ncol(M0),ncol(M0))
#         
#         # Initialize three matrices
#         Bold <- matrix(0,dim(U0)[1], dim(U0)[2])
#         Cold <- matrix(0,dim(U0)[1], dim(U0)[2])
#         etaold <- matrix(0,dim(U0)[1], dim(U0)[2])
#         
#         # The MAIN loop of admm method
#         step <- 0    
#         start_time <- Sys.time()
#         
#         repeat{
#           step <- step + 1
#           
#           # Update B
#           U <- U0 - etaold + gamma * Cold
#           out_B <- updateB_oracle(M, U, lambda1, opts)
#           Bnew <- out_B$Bnew
#           jerr <- out_B$jerr
#           if(jerr != 0) break
#           
#           # Update C
#           Cnew <- updateC_oracle(Bnew, lambda2, gamma, etaold)
#           
#           # Update U
#           etanew <- etaold + gamma * (Bnew - Cnew)
#           
#           # Code 1: success
#           if(max(abs(Bnew - Bold)) < eps){
#             jerr <- 1
#             break
#           }
#           # Code 404: then maximal iteration is reached
#           if(step > maxit){
#             jerr <- 404
#             if(warning.it >= 3) warning('Maximal iteration is reached.')
#             break
#           }
#           Bold <- Bnew
#           Cold <- Cnew
#           etaold <- etanew
#         }# End of repeat 
#         end_time <- Sys.time()  # The time for each repeat
#         time <- difftime(end_time, start_time, units = "secs")
#         # Code < -10000: non-sparse matrix
#         if(jerr < -10000){
#           if(warning.it >= 3) warning(paste0('The ',i, 'th lam1 is too small, no execution for the rest parameters lam2 and gamma.'))
#           break
#         }
#         # Code 1: success, save the matrix and the related information.
#         if(jerr==1){
#           index <- (i-1)*n2*n3 + (j-1)*n2 + k
#           nlam_cvg <- nlam_cvg + 1
#           B_l[[index]] <- Bnew
#           step_l[index] <- step
#           time_l[index] <- time
#           if(is.null(d)) rank <- rank_func(Cnew, thrd = 1e-3)
#           else rank <- d
#           rank_l[index] <- rank
#           # save the singular values of each candidates matrix B and C
#           sv_list_B[[index]] <- svd(Bnew)$d
#           sv_list_C[[index]] <- svd(Cnew)$d
#           # Cut and select the left singular vector of Bnew
#           if(rank == 0){
#             beta <- matrix(0, nrow(Bnew), ncol(Bnew))
#           }else{
#             tmp <- svd(Bnew)$u[,1:rank, drop = FALSE]
#             vec <- as.vector(tmp)
#             vec[abs(vec) < 1e-3] <- 0
#             beta <- matrix(vec, nrow(tmp), ncol(tmp))
#           }
#           beta_l[[index]] <- beta
#           var_ind <- apply(beta, 1, function(x){any(x!=0)})
#           s_l[index] <- sum(var_ind)
#         }
#       }# End of lambda2
#       if(jerr < -10000) break
#     }# End of gam
#   }# End of lambda1
#   return(list(beta = beta_l, B = B_l, rank = rank_l, s = s_l, step = step_l, time = time_l, nlam = nlam_cvg, sv_list_B = sv_list_B, sv_list_C = sv_list_C))
# }
# 
# updateB_oracle <- function(M, U, lambda1, opts){
#   U <- t(U)
#   fit <- .Fortran("msda", obj = opts$nlam, opts$H, opts$nvars, as.double(M), as.double(U), opts$pf, opts$dfmax, opts$pmax, opts$nlam, opts$flmin, lambda1, opts$eps_inner, opts$maxit_inner, opts$sml, opts$verbose, nalam = opts$nalam, theta = opts$theta, itheta = opts$itheta, ntheta = opts$ntheta, alam = opts$alam, npass = opts$npass, jerr = opts$jerr)
#   if(fit$jerr != 0){return(list(Bnew = NULL, jerr = fit$jerr))} # Code: non-zero, abnormal result.
#   outlist <- formatoutput(fit, opts$maxit_inner, opts$pmax, opts$nvars, opts$H)
#   Bnew <- as.matrix(outlist$theta[[1]])
#   list(Bnew = Bnew, jerr = fit$jerr)
# }
# 
# updateC_oracle <- function(Bnew, lambda2, gamma, etaold){
#   Btemp <- Bnew + 1/gamma * etaold
#   svd_B <- svd(Btemp)
#   lamtemp <- pmax(0, svd_B$d-lambda2/gamma)
#   Cnew <- svd_B$u %*% diag(lamtemp, nrow = length(lamtemp), ncol = length(lamtemp)) %*% t(svd_B$v)
#   Cnew
# }

eval_dc_oracle <- function(Beta, x, y, yclass = NULL, H = H, categorical = FALSE){
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  ydummy <- sapply(unique(yclass)[2:length(unique(yclass))], function(i){
    as.numeric(yclass == i)
  })
  if(!is.list(Beta)){Beta <- list(Beta)}
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      1 - dcor(x %*% mat, ydummy)
    }
  })
  return(result)
}
