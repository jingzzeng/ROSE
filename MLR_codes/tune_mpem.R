tune_mpem = function(X, Y, K, sig, fortran=TRUE, eta.init=NULL, seqlamb, iter.max=50, print=FALSE){
  n = dim(X)[1]
  p = dim(X)[2]
  minbic = 10^8
  
  for (ilamb in 1:length(seqlamb)){
    #print(seqlamb[ilamb])
    
    t1 = Sys.time()
    obj = mpem(X, Y, K, sig, fortran, eta.init, lambda=seqlamb[ilamb], iter.max=iter.max, print=print)
    obj_time = difftime(Sys.time(), t1, units="secs")
    
    beta = obj$beta
    sig_est = obj$sig
    #print(sig_est)
    
    
    
    # for (i in 1:n){
    #   tmp = 0
    #   logf1 = -(Y[i]-sum(X[i,]*beta[,1]))^2/(2*sig) - log(2*pi*sig)/2
    #   for (k in 1:K){
    #     logfkof1 = ((Y[i]-sum(X[i,]*beta[,1]))^2 - (Y[i]-sum(X[i,]*beta[,k]))^2)/(2*sig)
    #     fkoverf1 = exp(logfkof1)
    #     tmp = tmp + obj$wgt[k]*fkoverf1
    #   }
    #   crit = crit+log(tmp)+logf1
    # }
    
    loglk = loglk(X, Y, beta, obj$wgt, sig=sig_est)
    
    bic = -2*loglk + log(n)*sum(beta!=0)
    
    
    if (bic < minbic){
      minbic = bic
      opt_lamb = seqlamb[ilamb]
      opt_beta = obj$beta
      opt_id = obj$id
      opt_iter = obj$niter
      opt_time = obj_time
      opt_wgt = obj$wgt
      opt_sig = obj$sig
      opt_eta = obj$eta
    }
  }
  
  outlist = list(opt_bic=minbic,opt_lamb=opt_lamb,opt_beta=opt_beta,opt_wgt=opt_wgt,
                 opt_id=opt_id,opt_iter=opt_iter,opt_time=opt_time,sig=opt_sig,eta=opt_eta)
  return(outlist)
}


tune_mpem_K = function(X, Y, rangeK ,sig, fortran=TRUE, eta.init=NULL, seqlamb, iter.max=50, print=FALSE){
  n = dim(X)[1]
  p = dim(X)[2]
  bick=rep(0,length(rangeK))
  for (mk in 1: length(rangeK)) {
    K=rangeK[mk]
    initial = init_eta(X, Y, K, sig, beta.init="tensor", random.init=NULL)
    eta.init = initial$eta.est
    beta.init = initial$beta.est
    minbic = 10^8
    for (ilamb in 1:length(seqlamb)){
      #print(seqlamb[ilamb])
      
      t1 = Sys.time()
      obj = mpem(X, Y, K, sig, fortran, eta.init, lambda=seqlamb[ilamb], iter.max=iter.max, print=print)
      obj_time = difftime(Sys.time(), t1, units="secs")
      
      beta = obj$beta
      sig_est = obj$sig
      #print(sig_est)
      
      
      
      # for (i in 1:n){
      #   tmp = 0
      #   logf1 = -(Y[i]-sum(X[i,]*beta[,1]))^2/(2*sig) - log(2*pi*sig)/2
      #   for (k in 1:K){
      #     logfkof1 = ((Y[i]-sum(X[i,]*beta[,1]))^2 - (Y[i]-sum(X[i,]*beta[,k]))^2)/(2*sig)
      #     fkoverf1 = exp(logfkof1)
      #     tmp = tmp + obj$wgt[k]*fkoverf1
      #   }
      #   crit = crit+log(tmp)+logf1
      # }
      
      loglk1 = loglk(X, Y, beta, obj$wgt, sig=sig)
      
      bic = -2*loglk1 + log(n)*sum(beta!=0)
      
      
      if (bic < minbic){
        minbic = bic
        opt_lamb = seqlamb[ilamb]
        opt_beta = obj$beta
        opt_id = obj$id
        opt_iter = obj$niter
        opt_time = obj_time
        opt_wgt = obj$wgt
        opt_sig = obj$sig
        opt_eta = obj$eta
      }
      mpem_ind = which(opt_beta[,1]!=0)
      mpem3_fit = pem(X[,mpem_ind], Y, K, sig, group=FALSE, lambda=0, eta.init=opt_eta, iter.max=20)
      mpem3_beta = matrix(0,p,K)
      mpem3_beta[mpem_ind,] = mpem3_fit$beta
    }
    loglk1 = loglk(X, Y, mpem3_beta, obj$wgt, sig=sig)
    bic = -2*loglk1 + log(n)*sum(mpem3_beta!=0)
    bick[mk]=bic
  }
  opt_K=rangeK[which.min(bick)]
  outlist = list(opt_bic=minbic,opt_lamb=opt_lamb,opt_beta=opt_beta,opt_wgt=opt_wgt,
                 opt_id=opt_id,opt_iter=opt_iter,opt_time=opt_time,sig=opt_sig,eta=opt_eta,opt_K=opt_K,bic_K=bick)
  return(outlist)
}


tune_mpem_K1 = function(X, Y, rangeK ,sig, fortran=TRUE, eta.init=NULL, seqlamb, iter.max=50, print=FALSE){
  n = dim(X)[1]
  p = dim(X)[2]
  bick=rep(0,length(rangeK))
  for (mk in 1: length(rangeK)) {
    K=rangeK[mk]
    initial = init_eta(X, Y, K, sig, beta.init="tensor", random.init=NULL)
    eta.init = initial$eta.est
    beta.init = initial$beta.est
    minbic = 10^8
    for (ilamb in 1:length(seqlamb)){
      #print(seqlamb[ilamb])
      
      t1 = Sys.time()
      obj = mpem(X, Y, K, sig, fortran, eta.init, lambda=seqlamb[ilamb], iter.max=iter.max, print=print)
      obj_time = difftime(Sys.time(), t1, units="secs")
      
      beta = obj$beta
      sig_est = obj$sig
      #print(sig_est)
      
      
      
      # for (i in 1:n){
      #   tmp = 0
      #   logf1 = -(Y[i]-sum(X[i,]*beta[,1]))^2/(2*sig) - log(2*pi*sig)/2
      #   for (k in 1:K){
      #     logfkof1 = ((Y[i]-sum(X[i,]*beta[,1]))^2 - (Y[i]-sum(X[i,]*beta[,k]))^2)/(2*sig)
      #     fkoverf1 = exp(logfkof1)
      #     tmp = tmp + obj$wgt[k]*fkoverf1
      #   }
      #   crit = crit+log(tmp)+logf1
      # }
      
      loglk1 = loglk(X, Y, beta, obj$wgt, sig=sig_est)
      
      bic = -2*loglk1 + log(n)*sum(beta!=0)
      
      
      if (bic < minbic){
        minbic = bic
        opt_lamb = seqlamb[ilamb]
        opt_beta = obj$beta
        opt_id = obj$id
        opt_iter = obj$niter
        opt_time = obj_time
        opt_wgt = obj$wgt
        opt_sig = obj$sig
        opt_eta = obj$eta
      }
      # mpem_ind = which(opt_beta[,1]!=0)
      # mpem3_fit = pem(X[,mpem_ind], Y, K, sig, group=FALSE, lambda=0, eta.init=opt_eta, iter.max=20)
      # mpem3_beta = matrix(0,p,K)
      # mpem3_beta[mpem_ind,] = mpem3_fit$beta
    }
    # loglk1 = loglk(X, Y, mpem3_beta, obj$wgt, sig=sig)
    # bic = -2*loglk1 + log(n)*sum(mpem3_beta!=0)
    bick[mk]=minbic
  }
  opt_K=rangeK[which.min(bick)]
  outlist = list(opt_bic=minbic,opt_lamb=opt_lamb,opt_beta=opt_beta,opt_wgt=opt_wgt,
                 opt_id=opt_id,opt_iter=opt_iter,opt_time=opt_time,sig=opt_sig,eta=opt_eta,opt_K=opt_K,bic_K=bick)
  return(outlist)
}
