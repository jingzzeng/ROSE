# source("PEM_utilities.R")
# source("power_utilities.R")
# source("powerMethod.R")


init_eta = function(X, Y, K, sig = 1, beta.init=NULL, random.init=NULL){
  #beta.init: input a matrix for initialization, or choose "tensor", "hddc"
  n = dim(X)[1]
  p = dim(X)[2]
  if(is.null(sig)) sig <- 1
  
  if(is.null(beta.init)&is.null(random.init)){
    stop("Input initialized beta or which random initialization method to use!")
  }
  
  if(!is.null(beta.init)){
    
    if(is.matrix(beta.init)){
      eta.est = est_eta(X, Y, beta.init, wgt=rep(1,K)/K, sig)
      
    } else if(beta.init=="tensor"){
      p = dim(X)[2]
      cvlasso = cv.glmnet(X,Y,family="gaussian",intercept=FALSE)
      lambda1 = cvlasso$lambda.1se
      #lambda1 = cvlasso$lambda.min * 0.8
      fit = glmnet(X,Y,lambda=lambda1,family="gaussian",intercept=FALSE)
      b_ini = matrix(coef(fit))[-1]
      ind = (1:p)[b_ini>=0.01]
      #ind = 1:5
      X_cut = X[,ind,drop=F]
      print(ind)
      
      #-----------------------------
      # Second-moment
      #-----------------------------
      mat_est = powerMethod_mat(X_cut, Y, K, linear=TRUE, sparse=FALSE)
      beta_m2_basis = mat_est$beta_est
      beta_m2_prj = beta_m2_basis%*%t(beta_m2_basis)
      #beta_m2 = norm_match(beta_m2_basis,beta[ind,])$beta_new
      
      #--------------------------
      # Third-moment with whiten
      #--------------------------
      power_est = powerMethod(X_cut, Y, K, linear=TRUE, whiten=TRUE, M2=mat_est$T_moment, N_iter=20, L_iter=50)
      beta_tensor = power_est$beta_est
      #beta_tensor = norm_match(beta_tensor,beta[ind,])$beta_new
      #print(beta_tensor)
      print("power estimation done")
      
      #-----------------------------------------------------
      # Project X to the subspace and decompose third-moment 
      #-----------------------------------------------------
      # basis_est = powerMethod_mat(X_cut, Y, K, linear=TRUE, sparse=FALSE)$beta_est
      # X_red = X_cut%*%basis_est
      # 
      # mat_est = powerMethod_mat(X_red, Y, K, linear=TRUE, sparse=FALSE)
      # power_est = powerMethod(X_red, Y, K, linear=TRUE, whiten=TRUE, M2=mat_est$T_moment, N_iter=20, L_iter=50)
      # beta_m3x = basis_est%*%power_est$beta_est
      # print(power_est$beta_est )
      # print(beta_m3x)
      
      eta.est = est_eta(X_cut, Y, beta_tensor, wgt=rep(1,K)/K, sig)
      
      beta.est = matrix(0,p,K)
      beta.est[ind,] = beta_tensor
      
      
    } else if(beta.init=="hddc"){
      p = dim(X)[2]
      cvlasso = cv.glmnet(X,Y,family="gaussian",intercept=FALSE)
      lambda1 = cvlasso$lambda.min * 0.05
      fit = glmnet(X,Y,lambda=lambda1,family="gaussian",intercept=FALSE)
      b_ini = matrix(coef(fit))[-1]
      ind = (1:p)[b_ini>=0.01]
      X_cut = cbind(X[,ind,drop=F],Y)
      prms1 = hddc(X_cut, K=K, algo="EM", init='kmeans')
      eta.est = prms1$posterior
      # print(ind)
      
    } else {
      stop("beta initialization input incorrect!")
    }
    
  }else{
    if(random.init == "1"){
      #--------------------------
      #Random initialization 1
      #--------------------------
      beta.est = matrix(runif(p*K,-0.1,0.1),p,K)
      eta.est = est_eta(X, Y, beta.est, wgt=rep(1,K)/K, sig)
    }
    if(random.init == "2"){
      #--------------------------
      #Random initialization 2
      #--------------------------
      Z.ini = sample(K,n,replace=TRUE)
      eta.est = matrix(0.1,n,K)
      for (i in 1:n) {
        eta.est[i,Z.ini[i]] = 0.9
      }
      eta.est = eta.est/(0.9+0.1*(K-1))
    }
    
    if(random.init == "3"){
      #--------------------------
      #Random initialization 3
      #--------------------------
      eta.est = matrix(1/K,n,K)
    }
    
  }

  if(exists("beta.est")){
    return_list = list(eta.est=eta.est, beta.est=beta.est)
  } else {
    return_list = list(eta.est=eta.est, beta.est = NULL)
  }
  
  return(return_list)
}



init_eta1 = function(X, Y, K, sig, beta.init=NULL, random.init=NULL){
  #beta.init: input a matrix for initialization, or choose "tensor", "hddc"
  
  n = dim(X)[1]
  
  if(is.null(beta.init)&is.null(random.init)){
    stop("Input initialized beta or which random initialization method to use!")
  }
  
  if(is.null(beta.init)==FALSE){
    
    if(is.matrix(beta.init)){
      eta.est = est_eta(X, Y, beta.init, wgt=rep(1,K)/K, sig)
      
    } else if(beta.init=="hddc"){
      p = dim(X)[2]
      cvlasso = cv.glmnet(X,Y,family="gaussian",intercept=FALSE)
      lambda1 = cvlasso$lambda.min * 0.05
      fit = glmnet(X,Y,lambda=lambda1,family="gaussian",intercept=FALSE)
      b_ini = matrix(coef(fit))[-1]
      ind = (1:p)[b_ini>=0.01]
      X_cut = cbind(X[,ind,drop=F],Y)
      prms1 = hddc(X_cut, K=K, algo="EM", init='kmeans')
      eta.est = prms1$posterior
      print(ind)
      
    } else if(beta.init=="tensor"){
      p = dim(X)[2]
      cvlasso = cv.glmnet(X,Y,family="gaussian",intercept=FALSE)
      lambda1 = cvlasso$lambda.1se
      #lambda1 = cvlasso$lambda.min * 0.8
      fit = glmnet(X,Y,lambda=lambda1,family="gaussian",intercept=FALSE)
      b_ini = matrix(coef(fit))[-1]
      ind = (1:p)[b_ini>=0.01]
      ind=1:10
      #ind = 1:5
      X_cut = X[,ind,drop=F]
      print(ind)
      
      #-----------------------------
      # Second-moment
      #-----------------------------
      mat_est = powerMethod_mat(X_cut, Y, K, linear=TRUE, sparse=FALSE)
      beta_m2_basis = mat_est$beta_est
      beta_m2_prj = beta_m2_basis%*%t(beta_m2_basis)
      #beta_m2 = norm_match(beta_m2_basis,beta[ind,])$beta_new
      
      #--------------------------
      # Third-moment with whiten
      #--------------------------
      power_est = powerMethod(X_cut, Y, K, linear=TRUE, whiten=TRUE, M2=mat_est$T_moment, N_iter=20, L_iter=50)
      beta_tensor = power_est$beta_est
      #beta_tensor = norm_match(beta_tensor,beta[ind,])$beta_new
      #print(beta_tensor)
      print("power estimation done")
      
      #-----------------------------------------------------
      # Project X to the subspace and decompose third-moment 
      #-----------------------------------------------------
      # basis_est = powerMethod_mat(X_cut, Y, K, linear=TRUE, sparse=FALSE)$beta_est
      # X_red = X_cut%*%basis_est
      # 
      # mat_est = powerMethod_mat(X_red, Y, K, linear=TRUE, sparse=FALSE)
      # power_est = powerMethod(X_red, Y, K, linear=TRUE, whiten=TRUE, M2=mat_est$T_moment, N_iter=20, L_iter=50)
      # beta_m3x = basis_est%*%power_est$beta_est
      # print(power_est$beta_est )
      # print(beta_m3x)
      
      eta.est = est_eta(X_cut, Y, beta_tensor, wgt=rep(1,K)/K, sig)
      
      beta.est = matrix(0,p,K)
      beta.est[ind,] = beta_tensor
      
      
    } else {
      stop("beta initialization input incorrect!")
    }
    
  }
  
  
  
  
  if(is.null(random.init)==FALSE){
    if(random.init==1){
      #--------------------------
      #Random initialization 1
      #--------------------------
      Z.ini = sample(K,n,replace=TRUE)
      eta.est = matrix(0.1,n,K)
      for (i in 1:n) {
        eta.est[i,Z.ini[i]] = 0.9
      }
      eta.est = eta.est/(0.9+0.1*(K-1))
    }
    
    if(random.init==2){
      #--------------------------
      #Random initialization 2
      #--------------------------
      eta.est = matrix(1/K,n,K)
    }
    
    if(random.init==3){
      #--------------------------
      #Random initialization 3
      #--------------------------
      beta.init = matrix(runif(p*K,-0.1,0.1),p,K)
      eta.est = est_eta(X, Y, beta.init, wgt=rep(1,K)/K, sig)
    }
    
  }
  
  if(exists("beta.est")){
    return_list = list(eta.est=eta.est, beta.est=beta.est)
  } else {
    return_list = list(eta.est=eta.est)
  }
  
  return(return_list)
}



