est_eta = function(X, Y, beta, wgt, sig = 1){
  if(is.null(sig)) sig <- 1
  n = nrow(X)
  K = ncol(beta)
  
  eta = matrix(0,n,K)
  for(i in 1:n){
    for(j in 1:K){
      rp = 0
      for(k in 1:K){
        if(k==j){
          temp = 0
        } else {
          temp = sum((beta[,k]-beta[,j])*X[i,])*(Y[i]-sum((beta[,k]+beta[,j])*X[i,])/2)/sig
        }
        rp = rp + wgt[k]*exp(temp)
      }
      eta[i,j] = wgt[j] / rp
    }
  }
  eta[which(is.na(eta))] = 0
  
  return(eta)
}


pmse = function(Xnew,Ynew,beta,wgt,sig){
  
  n = nrow(Xnew)
  eta.new = est_eta(Xnew,Ynew,beta,wgt,sig)
  Z = apply(eta.new, 1, which.max)
  
  pmse = 0
  for (i in 1:n) {
    pmse = pmse + (Ynew[i]-sum(Xnew[i,]*beta[,Z[i]]))^2
  }
  pmse = pmse/n
  return(pmse)
}





loglk = function(X, Y, beta, wgt, sig){
  
  n = nrow(X)
  K = ncol(beta)
  
  crit = 0
  for (i in 1:n){
    tmp = 0
    for (k in 1:K){
      logfk = -(Y[i]-sum(X[i,]*beta[,k]))^2/(2*sig) - log(2*pi*sig)/2
      tmp = tmp + wgt[k]*exp(logfk)
    }
    crit = crit+log(tmp)
  }
  
  return(crit)
}



cluster_err = function(K, idx, id_est) {
  # K is number of clusters
  # id_est is the estimate label
  # idx is the true label
  # return error rate

  perK = combinat::permn(K)

  n = length(idx)

  K_pred = perK[[1]]
  id_pred = K_pred[id_est]
  for(t in 2:length(perK)) {
    K_now = perK[[t]]
    id_now = K_now[id_est]
    if(sum(id_now!=idx)<sum(id_pred!=idx)){
      K_pred = K_now
      id_pred = id_now
    }
  }

  id_err = sum(id_pred!=idx)/n*100

  return(list(cluster_err=id_err, K_pred=K_pred, id_pred=id_pred))
}



##############################################
##### Estimation error with permutation ######
##############################################
beta_err = function(beta, beta_est) {
  K = ncol(beta)
  beta_err = Inf
  
  perK = combinat::permn(K)
  for (i in 1:length(perK)) {
    beta_now = beta_est[,perK[[i]]]
    err_now = norm(beta-beta_now, type="F")
    if(err_now < beta_err){
      beta_err = err_now
      beta_perm = beta_now
    }
  }
  
  outlist = list(beta_err=beta_err, beta_perm=beta_perm)
  return(outlist)
}



####################################################
##### MM algorithm for beta estimation in PEM ######
####################################################
MLR=function(Sigma,rho,K,B,lambda,n.iter=100,eps=1e-3){
  p=dim(Sigma[[1]])[1]
  beta0=beta=beta1=matrix(B,p*K,1)
  U=rep(0,K*p)
  gamma=rep(0,K)
  for (k in 1:K){
    U[((k-1)*p+1):(k*p)]=-(Sigma[[k]]%*%beta[((k-1)*p+1):(k*p)]-rho[[k]])
  }
  U=matrix(U,K*p,1)
  #U=matrix(-c(Sigma1%*%beta[1:p]-rho1,Sigma2%*%beta[(1+p):(2*p)]-rho2),2*p,1)
  gamma=rep(0,K)
  for (iter in 1:n.iter){
    for (j in 1:p){
      for (k in 1:K){
        sigma=Sigma[[k]]
        gamma[k]=sigma[j,j]
      }
      r=1*max(gamma)+0.001
      ind=seq(0,K-1,1)*p
      ind=ind+j
      beta1[ind,]=1/r*(U[ind,]+r*beta[ind,])*ifelse(norm(matrix(U[ind,]+r*beta[ind,],K,1),'f')>lambda, 1-lambda/norm(matrix(U[ind,]+r*beta[ind,],K,1),'f'),0)
      for (k in 1:K){
        U[((k-1)*p+1):(k*p)]=U[((k-1)*p+1):(k*p)]+Sigma[[k]][,j]*(beta[(k-1)*p+j]- beta1[(k-1)*p+j] )
      }
      # U1=U[1:p]+Sigma1[,j]*(beta[j]-beta1[j])
      # U2=U[(p+1):(2*p)]+Sigma2[,j]*(beta[p+j]-beta1[p+j])
      # U=matrix(c(U1,U2),,1)
      beta[ind,]=beta1[ind,]
    }
    if (norm(beta-beta0,'f')<eps){break}
    beta0=beta
  }
  #b1=beta[1:p,];b2=beta[(p+1):(2*p),]
  #fvalue=t(b1)%*%Sigma1%*%b1+t(b2)%*%Sigma2%*%b2-2*t(b1)%*%rho1-2*t(b2)%*%rho2
  b=matrix(beta,p,K)
  list(beta=b)
}

perm.out <- function(K, beta, beta.est){
  perm.set <- rbind(1:K, allPerms(1:K))
  beta.error.pool <- apply(perm.set, 1, function(a){sqrt(sum((beta.est[,a,drop=FALSE] - beta)^2))})
  id.perm <- which.min(beta.error.pool)
  sele.perm <- perm.set[id.perm,]
  beta.error.perm <- beta.error.pool[id.perm]
  
  list(perm = sele.perm, beta.error = beta.error.perm)
}

