formatoutput <- function(fit, maxit, pmax, p, H, warning.it = 0) {
  nalam <- fit$nalam
  ntheta <- fit$ntheta[seq(nalam)]
  nthetamax <- max(ntheta)
  lam <- fit$alam[seq(nalam)]
  theta_vec <- fit$theta
  errmsg <- err(fit$jerr, maxit, pmax)  ### error messages from fortran
  if(warning.it >= 3) switch(paste(errmsg$n), `1` = stop(errmsg$msg, call. = FALSE), `-1` = cat(errmsg$msg))
  if(nthetamax > 0){
    ja <- fit$itheta[seq(nthetamax)]
    theta <- lapply(seq_len(nalam), function(i){
      tmp <- theta_vec[(pmax * H * (i-1) + 1):(pmax * H * i)]
      a <- matrix(tmp, pmax, H, byrow = TRUE)[seq(nthetamax), , drop = FALSE]
      theta_i <- matrix(0, p, H)
      theta_i[ja,] <- a
      theta_i
    })
  }
  else{
    theta <- lapply(seq(nalam), function(x){matrix(0, p, H)})
  }
  list(theta = theta, lambda = lam)
}

err <- function(n, maxit, pmax) {
  if (n == 0) 
    msg <- ""
  if (n > 0) {
    # fatal error
    if (n < 7777) 
      msg <- "Memory allocation error; contact package maintainer"
    if (n == 10000) 
      msg <- "All penalty factors are <= 0"
    n <- 1
    msg <- paste("in the fortran code -", msg)
  }
  if (n < 0) {
    # non fatal error
    if (n > -10000) 
      msg <- paste("Convergence for ", -n, "th lambda value not reached after maxit=", 
                   maxit, " iterations; solutions for larger lambdas returned.\n", 
                   sep = "")
    if (n < -10000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds pmax=", 
                   pmax, " at ", -n - 10000, "th lambda value; solutions for larger lambdas returned.\n", 
                   sep = "")
    if (n < -20000) 
      msg <- paste("Number of nonzero coefficients along the path exceeds dfmax=", 
                   pmax, " at ", -n - 20000, "th lambda value; solutions for larger lambdas returned.\n", 
                   sep = "")
    n <- -1
  }
  list(n = n, msg = msg)
}

lamfix <- function(lam){
  llam <- log(lam)
  if(length(llam) >= 3){lam[1] <- exp(2 * llam[2] - llam[3])}
  lam
}