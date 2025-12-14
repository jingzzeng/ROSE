ABESS_EM <- function(x, y, s = NULL, s.max = 15, K = NULL, K.max = 5, sigma = NULL, adaptive.sigma = TRUE, beta.init = NULL, eta.init = NULL, w.init = NULL, beta.init.type = c("tensor", "hddc"), random.init.type = c("1","2","3"), maxiter = 100, tol = 1e-3, tau = 1e-3){
  
  beta.init.type <- match.arg(beta.init.type)
  random.init.type <- match.arg(random.init.type)
  
  n <- NROW(x)
  p <- NCOL(x)
  
  # ------- Compute distance correlation for initializing the active set ------- #
  dist_cor <- sapply(seq_len(p), function(i){
    dcor(y, x[,i])
  })
  ord <- order(dist_cor, decreasing = TRUE)
  # ---------------------------------------------------------------------------- #
  
  # ---------- Initialization of w, beta, and eta ----------- #
  if(!is.null(K)){
    if((K != round(K)) || (K == 1)) stop("K should be some integer larger than 1.")
    K.list <- K
  }else{
    if((K.max != round(K.max)) || (K.max == 1)) stop("K.max should be some integer larger than 1.")
    K.list <- 2:K.max
  }
  if(!is.null(s)){
    s.list <- s
  }else{
    s.list <- 1:s.max
  }
  
  iter.tune <- rep(NA, length(K.list))
  sigma.list <- rep(NA, length(K.list))
  eta.list <- w.list <- pred.list <- vector("list", length = length(K.list))
  if(!is.null(beta.init) && !is.list(beta.init)) beta.init <- list(beta.init)
  if(!is.null(eta.init) && !is.list(eta.init)) eta.init <- list(eta.init)
  if(!is.null(w.init) && !is.list(w.init)) w.init <- list(w.init)
  
  if((length(K.list) == 1) && (length(s.list) == 1)){
    # If only one candidate for s and K is provided, no tuning is needed.
    K.opt <- K.list[1]
    s.opt <- s.list[1]
    eta.opt <- NULL
    w.opt <- NULL
    pred.opt <- NULL
    iter.opt <- NULL
  }else{
    ## If w.init is null, randomly generate it from Unif(0,1)
    if(is.null(w.init)){
      w.init <- vector("list", length = length(K.list))
      for(i in seq_along(K.list)){
        K.i <- K.list[i]
        w <- runif(K.i)
        w.init[[i]] <- w/sum(w)
      }
    }
    ## If eta.init is null, compute it from beta.init
    if(is.null(eta.init)){
      eta.init <- vector("list", length = length(K.list))
      ## If beta.init is null, use init_eta function to initialize eta and beta.
      if(is.null(beta.init)){
        for(i in seq_along(K.list)){
          K.i <- K.list[i]
          attempt.ini <- 1
          while (attempt.ini <= 10){
            try.out <- tryCatch({
              initial <- init_eta(x, y, K.i, sigma, beta.init=beta.init.type, random.init = random.init.type); NULL}, 
              error = function(e){message(message(sprintf("The %d-th attempt in initialization fails; try a new initialization.", attempt.ini)))
                message(paste0(e, "\n")); e})
            if(!inherits(try.out, "error")){
              break
            }
            attempt.ini <- attempt.ini + 1
          }
          if(attempt.ini > 10){
            stop("Multiple attempts of initialization fail. Need to consider another way of initialization.")
            return(NA)
          }
          if(!is.null(initial$beta.est)){beta.init[[i]] <- initial$beta.est}
          eta.init[[i]] <- initial$eta.est
        }
      }else{
        eta.init[[i]] <- est_eta(x, y, beta.init[[i]], w.init, sigma)
      }
    }
    if(is.null(beta.init)) beta.init <- vector("list", length = length(K.list))
    # ---------------------------------------------------------------------------- #
    
    
    # ------------ Main loop ------------- #
    s.old <- 0; s.opt <- s.max
    while(s.opt != s.old){ # If the optimal s is updated, do another round.
      s.old <- s.opt
      s.list <- 1:s.opt
      sic <- matrix(rep(NA, length(K.list) * length(s.list)), length(K.list), length(s.list))
      
      # ------------- Step 1: for each pair (K, s), run BESS_EM algorithm ----------- #
      for (i in seq_along(K.list)) {
        K.i <- K.list[i]
        
        # Run with the largest s
        j <- length(s.list)
        s.i <- s.list[j]
        cat(sprintf("(%d, %d)\n", K.i, s.i))
        
        # initialization
        A.init <- ord[1:(s.i)]
        A.init <- sort(A.init)
        w.i <- w.init[[i]]
        
        # Implement BESS_EM algorithm with given pair (K,s)
        attempt <- 1
        while (attempt <= 10){
          try.out <- tryCatch({
            beta.i <- beta.init[[i]]
            eta.i <- eta.init[[i]]
            output <- BESS_EM(x = x, y = y, A.init = A.init, K = K.i, sigma = sigma, adaptive.sigma = adaptive.sigma, w.init = w.i, beta.init = beta.i, eta.init = eta.i, maxiter = maxiter, tol = tol, tau = tau)
            NULL},
            error = function(e){
              message(sprintf("(K = %d, s = %d): The %d-th attempt in BESS_EM function fails; try another initial value.", K.i, s.i, attempt))
              message(paste0(e, "\n"))
              e
            })
          if(!inherits(try.out, "error")){
            break
          }else{
            attempt <- attempt + 1
            initial <- init_eta(x, y, K.i, sigma, beta.init=beta.init.type, random.init = random.init.type)
            if(!is.null(initial$beta.est)){beta.init[[i]] <- initial$beta.est}
            eta.init[[i]] <- initial$eta.est
          }
        }
        if(attempt > 10){
          warning("Multiple attempts of BESS_EM with different initialization fail. Need the manual investigation.")
          return(NA)
        }
        
        # record IC(s)
        L.s <- output$L
        A.s <- output$A
        sic[i,j] <- n * log(L.s) + 0.03 * length(A.s) * K.i * log(p) * log(log(n)) * (log(n))^2
        iter.tune[i] <- output$iter
        eta.list[[i]] <- eta <- output$eta
        w.list[[i]] <- output$w
        pred.list[[i]] <- output$pred
        sigma.list[i] <- output$sigma
        D <- vector("list", K.i)
        for(k in 1:K.i){
          D[[k]] <- diag(eta[,k])
        }
        
        # For the smaller s, implement BESS_EM algorithm with the fixed eta obtained from s = s_max.
        for (j in (length(s.list)-1):1) {
          s.i <- s.list[j]
          cat(sprintf("(%d, %d)\n", K.i, s.i))
          
          # initialization
          A.init <- ord[1:(s.i)]
          A.init <- sort(A.init)
          beta.i <- matrix(0, p, K.i)
          for(k in 1:K.i){
            # Initial beta using eta obtained from s = s_max
            beta.i[A.init, k] <- solve(t(x[,A.init,drop=FALSE]) %*% D[[k]] %*% x[,A.init,drop=FALSE]) %*% (t(x[,A.init,drop=FALSE]) %*% D[[k]] %*% y)
          }
          
          # Implement fast BESS_EM algorithm with fixed eta.
          output <- BESS_EM_fast(x = x, y = y, A.init = A.init, eta = eta, K = K.i, beta.init = beta.i, tau = tau)
          
          # record IC(s)
          L.s <- output$L
          A.s <- output$A
          sic[i,j] <- n * log(L.s) + 0.03 *  length(A.s) * K.i * log(p) * log(log(n))  * (log(n))^2
        }
      }
      
      # ------------- Step 2: select the optimal pair (K,s) ----------- #
      ind.opt <- which(sic == min(sic), arr.ind = TRUE)
      K.opt <- K.list[ind.opt[1]]
      s.opt <- s.list[ind.opt[2]]
      ## Record eta, w, pred, and iter obtained with the largest s in s.list since eta is fixed for smaller s.
      eta.opt <- eta.list[[ind.opt[1]]]
      w.opt <- w.list[[ind.opt[1]]]
      pred.opt <- pred.list[[ind.opt[1]]]
      sigma.opt <- sigma.list[ind.opt[1]]
      iter.opt <- iter.tune[ind.opt[1]]
      
      if(s.opt == 1) break
    }
  }
  
  # ------------- Step 3: rerun BESS_EM with s.opt  ----------- #
  if(is.null(eta.opt)){
    A.init <- ord[1:(s.opt)]
    A.init <- sort(A.init)
    w.i <- runif(K.opt)
    w.i <- w.i/sum(w.i)
    beta.i <- matrix(0, p, K.opt)
    beta.i[A.init,] <- matrix(rnorm(s.opt * K.opt), s.opt, K.opt)
    eta.i <- est_eta(x, y, beta.i, w.i, sigma)
    start.time <- Sys.time()
    output <- BESS_EM(x = x, y = y, A.init = A.init, K = K.opt, sigma = sigma, w.init = w.i, beta.init = beta.i, eta.init = eta.i, maxiter = maxiter, tol = tol, tau = tau)
    end.time <- Sys.time()
    time <- difftime(end.time, start.time, units = "secs")
    output <- c(list(K.opt = K.opt, s.opt = s.opt, time = time, eta.init = eta.init), output)
  }else{
    A.init <- ord[1:(s.opt)]
    A.init <- sort(A.init)
    D <- vector("list", K.opt)
    for(k in 1:K.opt){
      D[[k]] <- diag(eta.opt[,k])
    }
    start.time <- Sys.time()
    beta.i <- matrix(0, p, K.opt)
    for(k in 1:K.opt){
      beta.i[A.init, k] <- solve(t(x[,A.init,drop=FALSE]) %*% D[[k]] %*% x[,A.init,drop=FALSE]) %*% (t(x[,A.init,drop=FALSE]) %*% D[[k]] %*% y)
    }
    output <- BESS_EM_fast(x = x, y = y, A.init = A.init, eta = eta.opt, K = K.opt, beta.init = beta.i, tau = tau)
    end.time <- Sys.time()
    time <- difftime(end.time, start.time, units = "secs")
    output <- c(list(w = w.opt, pred = pred.opt, sigma = sigma.opt, iter = iter.opt, sic = sic, iter.tune = iter.tune, K.opt = K.opt, s.opt = s.opt, time = time, eta.init = eta.init), output)
  }
  
  output
}

BESS_EM_fast <- function(x, y, A.init, eta, K, beta.init = NULL, tau = 1e-3, h.max = length(A.init), beta.init.type = c("tensor", "hddc"), random.init.type = c("1","2","3")){
  A <- sort(A.init)
  A <- A.init
  D <- vector("list", K)
  for(k in 1:K){
    D[[k]] <- diag(eta[,k])
  }
  
  switch.out <- switch(x, y, K, eta, D, A, beta.init, h.max, tau)
  A <- switch.out$A
  L <- switch.out$L
  beta <- switch.out$beta
    
  out <- list(beta = beta, eta = eta, A = A, A.init = A.init, L = L)
  out  
}

BESS_EM <- function(x, y, A.init, K, sigma = NULL, adaptive.sigma = TRUE, w.init = NULL, beta.init = NULL, eta.init = NULL, maxiter = 20, tol = 1e-3, tau = 1e-3, h.max = length(A.init)){
  # n <- NROW(x)
  p <- NCOL(x)
  A.t <- sort(A.init)
  A.t <- A.init
  beta.t <- beta.init
  eta.t <- eta.init
  w.t <- w.init
  D.t <- vector("list", K)
  cvg.path <- c()
  
  if(!is.null(sigma)) adaptive.sigma <- FALSE
  if(isTRUE(adaptive.sigma)) sigma <- 1
  # Initialization
  if(is.null(eta.t)){
    if(is.null(beta.t)){
      initial <- init_eta(x, y, K, sigma, beta.init=beta.init.type, random.init = random.init.type)
      beta.t <- initial$beta.est
      eta.t <- initial$eta.est
    }else{
      eta.t <- est_eta(x, y, beta.t, w.t, sigma)
    }
  }
  for(k in 1:K){
    D.t[[k]] <- diag(eta.t[,k])
  }
  if(is.null(beta.t)) beta.t <- matrix(rnorm(p*K), p, K)
  
  ##
  sigma.path <- c(sigma)
  ##
  # ------------------ loop: EM algorithm ----------------- #
  for (t in 1:maxiter){
    
    # step 1: update w
    w.t <- colMeans(eta.t)
    
    # step 2: update beta.tilde.
    beta.tilde <- matrix(0, p, K)
    for(k in 1:K){
      tryCatch({
        beta.tilde[A.t, k] <- solve(t(x[,A.t,drop=FALSE]) %*% D.t[[k]] %*% x[,A.t,drop=FALSE]) %*% (t(x[,A.t,drop=FALSE]) %*% D.t[[k]] %*% y)
      },
      error = function(e){
        message(sprintf("Error in the update of the %d-th beta.tilde", k))
        message(paste0(e, "\n"))
      })
    }
    
    # step 3: switch variables and update A set
    switch.out <- switch(x, y, K, eta.t, D.t, A.t, beta.tilde, h.max, tau)
    
    A.t <- switch.out$A
    L.t <- switch.out$L
    beta.tmp <- switch.out$beta
    
    # step 4: check the convergence for EM algorithm.
    cvg <- sqrt(sum((beta.tmp - beta.t)^2))
    cvg.path <- c(cvg.path, cvg)
    beta.t <- beta.tmp # update beta
    if(cvg < tol){
      break
    }
    
    # (extra step: update sigma)
    if(isTRUE(adaptive.sigma)) sigma <- update.sigma(x, y, beta.t, eta.t)
    
    ##
    sigma.path <- c(sigma.path, sigma) # Take a look as the convergence path of sigma.
    ##
    
    # step 5: update eta
    eta.t <- est_eta(x, y, beta.t, w.t, sigma)
    for(k in 1:K){
      D.t[[k]] <- diag(eta.t[,k])
    }
  } 
  # ------------------ end of loop ----------------- #
  
  pred <- apply(eta.t, 1, which.max) # prediction for each sample
  out <- list(w = w.t, beta = beta.t, eta = eta.t, A = A.t, A.init = A.init, L = L.t, pred = pred, iter = t, cvg.path = cvg.path, sigma = sigma, sigma.path = sigma.path)
  out  
}

switch <- function(x, y, K, eta, D, A, beta, h.max, tau){
  n <- NROW(x)
  p <- NCOL(x)
  # initialization
  I <- setdiff(1:p, A)
  L <- sum((y - x[,A,drop=FALSE] %*% beta[A,,drop=FALSE])^2 * eta)/(2*n)
  if(identical(A, 1:p)) return(list(A = A, I = I, beta = beta, L = L))
  
  L.old <- L
  x.sq <- x^2
  betasq.slice <- (beta[A,,drop=FALSE])^2
  xi <- ((t(x.sq[,A,drop=FALSE]) %*% eta) * betasq.slice) %*% rep(1, K) /(2*n)
  
  eps <- y %*% t(rep(1,K)) - x[,A,drop=FALSE] %*% beta[A,,drop=FALSE]
  d <- t(x[,I,drop=FALSE]) %*% (eta * eps)
  d.sq <- d^2
  denom <- t(x.sq[,I,drop=FALSE]) %*% eta
  zeta <- (d.sq/denom) %*% rep(1,K) / (2*n)
  
  for (h in 1:h.max){
    sele.out <- sele(x, y, x.sq, eta, K, A, I, xi, zeta, h)
    A.new <- sele.out$A
    I.new <- sele.out$I
    beta.new <- matrix(0, p, K)
    for(k in 1:K){
      beta.new[A.new, k] <- solve(t(x[,A.new,drop=FALSE]) %*% D[[k]] %*% x[,A.new,drop=FALSE]) %*% (t(x[,A.new,drop=FALSE]) %*% D[[k]] %*% y)
    }
    L.new <- sum((y - x[,A.new,drop=FALSE] %*% beta.new[A.new,,drop=FALSE])^2 * eta)/(2*n)
    
    # if loss decreases, then switch
    if(L.new < L.old){
      A.old <- A.new
      I.old <- I.new
      beta.old <- beta.new
      L.old <- L.new
    }
  }
  
  # if loss decreases by at least tau, then update
  if(L - L.old > tau){
    A <- A.old
    I <- I.old
    beta <- beta.old
    L <- L.old
  }
  
  out <- list(A = A, I = I, beta = beta, L = L)
  out
}

sele <- function(x, y, x.sq, eta, K, A, I, xi, zeta, h){
  p <- NCOL(x)
  
  rm.set <- A[order(xi)[1:h]]
  add.set <- I[order(zeta, decreasing = TRUE)[1:h]]
  
  A.new <- union(setdiff(A, rm.set), add.set)
  A.new <- sort(A.new)
  I.new <- setdiff(1:p, A.new)
  
  list(A = A.new, I = I.new, xi = xi, zeta = zeta)
}

# ABESS_EM2 <- function(x, y, s = NULL, s.max = 15, K = NULL, K.max = 5, sigma, maxiter = 100, tol = 1e-3, tau = 1e-3){
#   
#   n <- NROW(x)
#   p <- NCOL(x)
#   
#   # ------------- compute distance correlation ----------- #
#   dist_cor <- sapply(seq_len(p), function(i){
#     dcor(y, x[,i])
#   })
#   ord <- order(dist_cor, decreasing = TRUE)
#   
#   # ------------- Step 1: for each pair (K, s), run BESS_EM algorithm ----------- #
#   if(!is.null(K)){
#     if((K != round(K)) || (K == 1)) stop("K should be some integer larger than 1.")
#     K.list <- K
#   }else{
#     if((K.max != round(K.max)) || (K.max == 1)) stop("K.max should be some integer larger than 1.")
#     K.list <- 2:K.max
#   }
#   if(!is.null(s)){
#     s.list <- s
#   }else{
#     s.list <- 1:s.max
#   }
#   
#   # If there is more than one pair of (K,s), we need to tune the parameters.
#   # iter.tune <- rep(NA, length(K.list))
#   # eta.list <- vector("list", length = length(K.list))
#   # w.list <- vector("list", length = length(K.list))
#   # pred.list <- vector("list", length = length(K.list))
#   if((length(K.list) == 1) && (length(s.list) == 1)){
#     K.opt <- K.list[1]
#     s.opt <- s.list[1]
#     eta.opt <- NULL
#     w.opt <- NULL
#     pred.opt <- NULL
#     iter.opt <- NULL
#   }else{
#     # s.old <- s.max; s.opt <- 0
#     # while(s.opt != s.old){
#       # s.old <- s.max
#       sic <- matrix(rep(NA, length(K.list) * length(s.list)), length(K.list), length(s.list))
#       for (i in seq_along(K.list)) {
#         K.i <- K.list[i]
#         
#         for (j in seq_along(s.list)) {
#           # j <- length(s.list)
#           s.i <- s.list[j]
#           cat(sprintf("(%d, %d)\n", K.i, s.i))
#           
#           # initialization
#           A.init <- ord[1:(s.i)]
#           A.init <- sort(A.init)
#           w.init <- runif(K.i)
#           w.init <- w.init/sum(w.init)
#           beta.init <- matrix(0, p, K.i)
#           beta.init[A.init,] <- matrix(rnorm(s.i*K.i), s.i, K.i)
#           
#           # Implement BESS_EM algorithm with given pair (K,s)
#           start.time <- Sys.time()
#           output <- BESS_EM(x = x, y = y, A.init = A.init, K = K.i, sigma = sigma, w.init = w.init, beta.init = beta.init, maxiter = maxiter, tol = tol, tau = tau)
#           end.time <- Sys.time()
#           time <- difftime(end.time, start.time, units = "secs")
#           
#           # record IC(s)
#           L.s <- output$L
#           A.s <- output$A
#           sic[i,j] <- n * log(L.s) + length(A.s) * K.i * log(p) * log(log(n))
#           # iter.tune[i] <- output$iter
#           # eta.list[[i]] <- eta <- output$eta
#           # w.list[[i]] <- output$w
#           # pred.list[[i]] <- output$pred
#           # D <- vector("list", K.i)
#           # for(k in 1:K.i){
#           #   D[[k]] <- diag(eta[,k])
#           # }
#         }
#         
#         # for (j in (length(s.list)-1):1) {
#         #   s.i <- s.list[j]
#         #   cat(sprintf("(%d, %d)\n", K.i, s.i))
#         #   
#         #   # initialization
#         #   A.init <- ord[1:(s.i)]
#         #   A.init <- sort(A.init)
#         #   beta.init <- matrix(0, p, K.i)
#         #   for(k in 1:K.i){
#         #     beta.init[A.init, k] <- solve(t(x[,A.init,drop=FALSE]) %*% D[[k]] %*% x[,A.init,drop=FALSE]) %*% (t(x[,A.init,drop=FALSE]) %*% D[[k]] %*% y)
#         #   }
#         #   
#         #   # Implement BESS_EM algorithm with given pair (K,s)
#         #   output <- BESS_EM_fast(x = x, y = y, A.init = A.init, eta = eta, K = K.i, beta.init = beta.init, tau = tau)
#         #   
#         #   # record IC(s)
#         #   L.s <- output$L
#         #   A.s <- output$A
#         #   sic[i,j] <- n * log(L.s) + length(A.s) * K.i * log(p) * log(log(n))
#         # }
#       }
#       
#       # ------------- Step 2: select the optimal pair (K,s) ----------- #
#       ind.opt <- which(sic == min(sic), arr.ind = TRUE)
#       K.opt <- K.list[ind.opt[1]]
#       s.opt <- s.list[ind.opt[2]]
#       # eta.opt <- eta.list[[ind.opt[1]]]
#       # w.opt <- w.list[[ind.opt[1]]]
#       # pred.opt <- pred.list[[ind.opt[1]]]
#       # iter.opt <- iter.tune[ind.opt[1]]
#       
#       # if(s.opt == 1) break
#       # 
#       # s.max <- s.opt
#       # s.list <- 1:s.max
#     # }
#   }
#   
#   # ------------- Step 3: rerun BESS_EM with s.opt  ----------- #
#   # if((is.null(eta.opt) || is.null(w.opt)) || is.null(pred.opt)){
#   A.init <- ord[1:(s.opt)]
#   A.init <- sort(A.init)
#   w.init <- runif(K.opt)
#   w.init <- w.init/sum(w.init)
#   beta.init <- matrix(0, p, K.opt)
#   beta.init[A.init,] <- matrix(rnorm(s.opt * K.opt), s.opt, K.opt)
#   start.time <- Sys.time()
#   output <- BESS_EM(x = x, y = y, A.init = A.init, K = K.opt, sigma = sigma, w.init = w.init, beta.init = beta.init, maxiter = maxiter, tol = tol, tau = tau)
#   end.time <- Sys.time()
#   time <- difftime(end.time, start.time, units = "secs")
#   output <- c(list(sic = sic, K.opt = K.opt, s.opt = s.opt, time = time), output)
#   # }else{
#   #   A.init <- ord[1:(s.opt)]
#   #   A.init <- sort(A.init)
#   #   D <- vector("list", K.opt)
#   #   for(k in 1:K.opt){
#   #     D[[k]] <- diag(eta.opt[,k])
#   #   }
#   #   beta.init <- matrix(0, p, K.opt)
#   #   for(k in 1:K.opt){
#   #     beta.init[A.init, k] <- solve(t(x[,A.init,drop=FALSE]) %*% D[[k]] %*% x[,A.init,drop=FALSE]) %*% (t(x[,A.init,drop=FALSE]) %*% D[[k]] %*% y)
#   #   }
#   #   output <- BESS_EM_fast(x = x, y = y, A.init = A.init, eta = eta.opt, K = K.opt, beta.init = beta.init, tau = tau)
#   #   output <- c(list(w = w.opt, pred = pred.opt, iter = iter.opt, sic = sic, iter.tune = iter.tune, K.opt = K.opt, s.opt = s.opt, time = time), output)
#   # }
#   
#   output
# }