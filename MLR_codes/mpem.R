#dyn.load("smlr.so")
# source("PEM_utilities.R")
#source("msda.R")

mpem = function(X, Y, K, sig=NULL, fortran=TRUE, eta.init=NULL, lambda=NULL, iter.max=50, stop=1e-3, print=FALSE){
  #X, n by p maxtrix of observed covariates
  #Y, a length n vector of univariate responses
  #K, number of heterogeneous groups
  #trueZ, a length n vector of true group labels
  #iter.extra, number of extra iterations to run if iter.max is achieved in order to select estimation with largest likelihood
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  if(is.null(sig)){
    sigmark = 1
  } else { sigmark = 0 }
  
  
  ## center data ##
  #X = scale(X,center=TRUE,scale=FALSE)
  
  #################################
  ###### Prepare for Fortran ######
  #################################
  nvars = as.integer(p)
  nobs = as.integer(n)
  nk = as.integer(K)
  dfmax = nobs
  pmax = nvars
  pf = rep(1,nvars)
  eps = 1e-04
  maxit = 1e+03
  sml = 1e-06
  ulam = lambda
  nlam = as.integer(1)

  pf = as.double(pf)
  if (length(pf) != nvars) 
    stop("The size of penalty factor must be same as the number of input variables")
  flmin = as.double(1)
  eps = as.double(eps)
  maxit = as.integer(maxit)
  verbose = as.integer(FALSE)
  sml = as.double(sml)
  
  
  ##########################
  ###### Penalized EM ######
  ##########################
  ests = list()
  
  
  #-----------------------
  # initialization eta
  #-----------------------
  if(is.null(eta.init)){
    stop("Input initialization for eta!")
  } else {
    ests$eta = eta.init
    ests$pi = colMeans(ests$eta)
  }
  
  sigma = t(X)%*%X/n
  #sigma = t(X)%*%X/n + 0.01*diag(p)
  #sigma = cov(X)
  #sigma = ARM(p,0.3)
  
  
  iter = 1
  ext.iter = 0
  while(iter<=iter.max){
    #print(iter)
    t0 = Sys.time()
    
    #lambda = lambda * 0.85
    
    #-------------------
    # M-step
    #-------------------
    rho.est = t(t(ests$eta)%*%(X*Y)/colSums(ests$eta))
    #rho.est = t(t(eta.est)%*%(X*Y)/n)
    #print(eta.est)
    #print(rho.est)
    
    # sigma = matrix(0,p,p)
    # Sigma = array(list(),K)
    # for (k in 1:K){
    #   Sigma[[k]] = matrix(0,p,p)
    #   for (i in 1:n){
    #     Sigma[[k]] = Sigma[[k]] + eta.est[i,k]*X[i,]%*%t(X[i,])/n
    #   }
    #   sigma = sigma + pi.est[k]*Sigma[[k]]
    # }
    
    if(fortran){
      delta = t(rho.est)
      
      fit = .Fortran("msda", obj=double(nlam), nk, nvars, as.double(sigma), as.double(delta), 
                     pf, dfmax, pmax, nlam, flmin, ulam, eps, maxit, sml, verbose, nalam=integer(1), 
                     theta=double(pmax*nk*nlam), itheta=integer(pmax), ntheta=integer(nlam), 
                     alam=double(nlam), npass=integer(1), jerr=integer(1))
      ests$beta = matrix(fit$theta, p, K, byrow=TRUE)
      
      #print(sum(beta_new[,1]!=0))
    } else {
      
      # sigma_true = diag(p)
      # for (i in 1:p){
      #   for (j in 1:p){
      #     sigma_true[i,j] = 0.3^abs(i-j)
      #   }
      # }
      
      beta_new = msda(sigma, rho.est, lambda)$beta
      
    }
    #print(beta_new)
    
    
    if(sigmark==1){
      epsilon = Y - X%*%ests$beta
      
      sighat = sum(ests$eta * epsilon^2)/n
      ests$sig = sighat
      #print(sig)
    } else {
      ests$sig = sig
    }
    
    
    #-------------------------------
    # E-step: calculate eta.est
    #-------------------------------
    ests$eta = est_eta(X, Y, ests$beta, ests$pi, ests$sig)
    #print(eta.est)
    ests$loglk = loglk(X, Y, ests$beta, ests$pi, ests$sig)
    #ests$loglk_pl = ests$loglk - lambda*sum(apply(ests$beta,1,norm,"2"))*n/(2*ests$sig)
    ests$pi = colMeans(ests$eta)
    #print(pi.est)
    
    if(iter>1){
      beta_err = norm(ests$beta-ests_old$beta,type="F") / norm(ests_old$beta,type="F")
      #print(norm(beta_new-beta_old,type="F"))
      #print(norm(beta_old,type="F"))
      if(is.na(beta_err)) beta_err = 0
      #beta_err = norm(beta_new-beta_old,type="F")
    } else {
      beta_err = Inf
    }
    
    ests$id = apply(ests$eta, 1, which.max)
    
    t_iter = difftime(Sys.time(), t0, units="secs")
    
    
    if (print) {
        cat(paste0("iter",iter),beta_err,ests$loglk,ests$sig,"\n",sep=" ")
    }


    
    
    if(iter>2){
      #if(ests_before$loglk>=(ests_old$loglk-1e-4) & ests_before$loglk>=(ests$loglk-1e-4) & (ests$loglk!=-Inf))
      #if(ests_before$loglk>=(ests_old$loglk+1) & ests_before$loglk>=(ests$loglk+1) & (ests$loglk!=-Inf))  
      if(ests_before$loglk>=(ests_old$loglk-1e-4) & ests_old$loglk>=(ests$loglk-1e-4) & (ests$loglk!=-Inf))
      {
        ests = ests_before
        break
      }
    }
    

    
    
    if(iter>1){
      ### record estimations including likelihood from 2 iterations before
      ests_before = ests_old
    }
    
    ests_old = ests
    
    
    iter = iter + 1   
    
      
  }
  
  
  outlist = list(wgt=ests$pi, eta=ests$eta,
                 beta=ests$beta, beta_err=beta_err,
                 id=ests$id, niter=iter-1, sig=ests$sig)
  
  return(outlist)
}
















