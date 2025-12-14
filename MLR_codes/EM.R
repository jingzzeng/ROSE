EM <- function(x, y, K, sigma = NULL, adaptive.sigma = TRUE, beta.init = NULL, eta.init = NULL, w.init = NULL, maxiter = 20, tol = 1e-3){
  p <- NCOL(x)
  n <- NROW(x)
  cvg.path <- c()
  if(!is.null(sigma)) adaptive.sigma <- FALSE
  if(isTRUE(adaptive.sigma)) sigma <- 1
  
  if(is.null(w.init)){
    w.t <- runif(K)
    w.t <- w.t/sum(w.t)
  }else{
    w.t <- w.init
  }
  if(is.null(beta.init)){
    beta.t <- matrix(rnorm(p * K), p, K)
  }else{
    beta.t <- beta.init
  }
  if(is.null(eta.init)){
    eta.t <- est_eta(x, y, beta.t, w.t, sig = 1)
  }else{
    eta.t <- eta.init
  }
  D.t <- vector("list", K)
  for(k in 1:K){
    D.t[[k]] <- diag(eta.t[,k])
  }
  
  # ------------------ loop: EM algorithm ----------------- #
  for (t in 1:maxiter){
    
    # step 1: update w
    w.t <- colMeans(eta.t)

    # step 2: update beta
    beta.tmp <- matrix(0, p, K)
    for(k in 1:K){
      beta.tmp[, k] <- solve(t(x) %*% D.t[[k]] %*% x) %*% (t(x) %*% D.t[[k]] %*% y)
    }

    # step 3: check the convergence for EM algorithm.
    cvg <- sqrt(sum((beta.tmp - beta.t)^2))
    cvg.path <- c(cvg.path, cvg)
    beta.t <- beta.tmp # update beta
    if(cvg < tol){
      break
    }
    
    # (extra step: update sigma)
    if(isTRUE(adaptive.sigma)) sigma <- update.sigma(x, y, beta.t, eta.t)
    
    # step 4: update eta
    for(k in 1:K){
      beta.diff <- beta.t[,-k,drop=FALSE] - beta.t[,k,drop=FALSE] %*% t(rep(1, K-1))
      beta.ave <- (beta.t[,-k,drop=FALSE] + beta.t[,k,drop=FALSE] %*% t(rep(1, K-1)))/2
      tmp <- exp((x %*% beta.diff) * (y - x %*% beta.ave) / sigma^2)  %*% w.t[-k]
      eta.t[,k] <- w.t[k]/(tmp + w.t[k])
      D.t[[k]] <- diag(eta.t[,k])
    }
  }
  # ------------------ end of loop ----------------- #

  L <- sum((y - x %*% beta.t)^2 * eta.t)/n
  beta <- beta.t
  pred <- apply(eta.t, 1, which.max) # prediction for each sample
  out <- list(w = w.t, beta = beta, eta = eta.t, pred = pred, iter = t, cvg.path = cvg.path, L = L, sigma = sigma)
  out
}