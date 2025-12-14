# source("PEM_utilities.R")

pem = function(X, Y, K, sig=NULL, group=TRUE, eta.init=NULL, lambda=NULL, iter.max=100, stop=1e-3, trueZ=NULL, print=FALSE){
  #X, n by p maxtrix of observed covariates
  #Y, a length n vector of univariate responses
  #K, number of heterogeneous groups
  #group, impose group sparsity, if FALSE, use general sparsity (lasso)
  #trueZ, a length n vector of true group labels
  
  n = dim(X)[1]
  p = dim(X)[2]
  
  if(is.null(sig)){
    sigmark = 1
  } else { sigmark = 0 }
  
  ## center data ##
  #X = scale(X,center=TRUE,scale=FALSE)
  
  ##########################
  ###### Penalized EM ######
  ##########################
  
  #-----------------------
  # initialization eta
  #-----------------------
  if(is.null(eta.init)){
    stop("Input initialization for eta!")
  } else {
    eta.est = eta.init
    pi.est = colMeans(eta.est)
  }
    
  #sigma = t(X)%*%X/n

  
  for(iter in 1:iter.max){
    t0 = Sys.time()
    
    #-------------------
    # M-step
    #-------------------
    beta_new = matrix(0,p,K)
    
    if(group){
      eta = eta.est
      rho = Sigma = vector('list',K)
      
      # rho1=rho2=matrix(0,p,1)
      # Sigma1=Sigma2=matrix(0,p,p)
      
      for (k in 1:K){
        rho[[k]] = matrix(0,p,1)
        Sigma[[k]] = matrix(0,p,p)
        for (i in 1:n){
          rho[[k]] = rho[[k]] + as.numeric(Y[i])*eta[i,k]*X[i,]/n
          Sigma[[k]] = Sigma[[k]] + eta[i,k]*X[i,]%*%t(X[i,])/n
        }
      }
      
      #obj1 = cv.MLR(X,Y,eta,lambda.max = 1,rate=0.02)
      obj = MLR(Sigma,rho,K,B=beta_new,lambda=lambda,n.iter=30)
      beta_new = obj$beta
    } else {
      for (k in 1:K) {
        cons = sqrt(eta.est[,k])
        Xnew = cons*X
        Ynew = cons*Y
        #cvlasso = cv.glmnet(Xnew,Ynew,family="gaussian",intercept=FALSE)
        #lambda1 = cvlasso$lambda.min
        fit = glmnet(Xnew,Ynew,lambda=lambda,family="gaussian",intercept=FALSE)
        
        beta_new[,k] = matrix(coef(fit))[-1]
        
      }
    }
    #print(beta_new)
    
    if(sigmark==1){
      epsilon = Y - X%*%beta_new
      
      sighat = sum(eta.est * epsilon^2)/n
      sig = sighat
      #print(sig)
    }

    
    
    #-------------------------------
    # E-step: calculate eta.est
    #-------------------------------
    eta.est = est_eta(X, Y, beta_new, pi.est, sig)
    pi.est = colMeans(eta.est)
    
    
    if(iter>1){
      beta_err = norm(beta_new-beta_old,type="F") / norm(beta_old,type="F")
      if(is.na(beta_err)) beta_err = 0
    } else {
      beta_err = Inf
    }
    
    id_fit = apply(eta.est, 1, which.max)
    
    t_iter = difftime(Sys.time(), t0, units="secs")
    
    
    if (print) {
      if(is.null(trueZ)){
        cat(paste0("iter",iter),beta_err,"\n",sep=" ")
      } else {
        id_err = cluster_err(K,trueZ,id_fit)$cluster_err
        cat(paste0("iter",iter),beta_err,id_err,"\n",sep=" ")
      }
    }
    
    if(beta_err<stop){
      break
    } else {
      beta_old = beta_new
    }
  }

  outlist = list(wgt=pi.est, eta=eta.est,
                 beta=beta_new, beta_err=beta_err,
                 id=id_fit, niter=iter-1, sig=sig)
  
  return(outlist)
}











