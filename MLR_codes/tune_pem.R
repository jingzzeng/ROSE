tune_pem = function(X, Y, K, sig, group=TRUE, eta.init=NULL, seqlamb, iter.max=20){
  n = dim(X)[1]
  p = dim(X)[2]
  minbic = 10^8
  
  for (ilamb in 1:length(seqlamb)){
    #print(seqlamb[ilamb])
    
    t1 = Sys.time()
    obj = pem(X, Y, K, sig, group, eta.init, lambda=seqlamb[ilamb], iter.max=iter.max)
    obj_time = difftime(Sys.time(), t1, units="secs")
    
    beta = obj$beta
    sig_est = obj$sig
    
    crit = 0
    
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
    
    for (i in 1:n){
      tmp = 0
      for (k in 1:K){
        logfk = -(Y[i]-sum(X[i,]*beta[,k]))^2/(2*sig_est) - log(2*pi*sig_est)/2
        tmp = tmp + obj$wgt[k]*exp(logfk)
      }
      crit = crit+log(tmp)
    }
    
    bic = -2*crit + log(n)*sum(beta!=0)
    
    if (bic < minbic){
      minbic = bic
      opt_lamb = seqlamb[ilamb]
      opt_beta = obj$beta
      opt_id = obj$id
      opt_iter = obj$niter
      opt_time = obj_time
      opt_sig = obj$sig
      ##
      opt_wgt = obj$wgt
    }
  }
  
  outlist = list(opt_bic=minbic,opt_lamb=opt_lamb,opt_beta=opt_beta,
                 opt_id=opt_id,opt_iter=opt_iter,opt_time=opt_time,sig=opt_sig,
                 opt_wgt = opt_wgt)
  
  return(outlist)
}



