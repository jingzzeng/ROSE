# source("power_utilities.R")

# powerMethod = function(X, y, K, GaussianInput=TRUE, linear=FALSE, whiten=TRUE, N_iter=100, L_iter=20, nu=0.1){
#   # X is observed predictors
#   # y is observed responses
#   # K is number of clusters
#   # If assume coffecient weights are orthogonal, set whiten=FALSE
#   
#   p = dim(X)[2]
#   beta_est = matrix(0,p,K)
#   
#   T_moment = moment3rd(y,X,linear=linear)
#   
#   if(whiten){
#     whitening = whiten(T_moment,r=K)
#     T_now = whitening$T_whitened
#     W = whitening$W
#     V = whitening$V
#     v = whitening$v
#   }
#   else{
#     T_now = T_moment 
#   }
# 
#   for (k in 1:K) {
#     
#     decomp = powerDecomp_rank1(T_now, N=N_iter, L=L_iter, nu=nu)
#     beta_temp = decomp$theta
#     
#     T_now = T_now - decomp$lamb * outer3d(beta_temp,beta_temp,beta_temp)
#     
#     if(whiten){
#       beta_est[,k] = MASS::ginv(t(W))%*%beta_temp / sqrt(v[k])
#     }
#     else{
#       beta_est[,k] = beta_temp
#     }
#     
#   }
#   
#   return_list = list(beta_est=beta_est)
#   return(return_list)
# }



powerMethod = function(X, y, K, GaussianInput=TRUE, linear=FALSE, whiten=TRUE, M2=NULL, N_iter=100, L_iter=20){
  # X is observed predictors
  # y is observed responses
  # K is number of clusters
  # If assume coffecient weights are orthogonal, set whiten=FALSE
  
  p = dim(X)[2]
  beta_est = matrix(0,p,K)
  
  T_moment = moment3rd(y,X,linear=linear)
  
  if(whiten){
    whitening = whiten(T_moment,r=K,V=M2)
    T_now = whitening$T_whitened
    W = whitening$W
    V = whitening$V
    v = whitening$v
  }
  else{
    T_now = T_moment 
  }
  
  
  if((whiten&linear)){
    pi.power = rep(0,K)
  }
  
  
  for (k in 1:K) {
    
    decomp = powerDecomp_rank1(T_now, N=N_iter, L=L_iter)
    beta_temp = decomp$theta
    lamb_temp = decomp$lamb
    
    T_now = T_now - lamb_temp * outer3d(beta_temp,beta_temp,beta_temp)
    
    if((whiten&linear)==TRUE){
      beta_est[,k] = MASS::ginv(t(W))%*%beta_temp * lamb_temp/3
      pi.power[k] = 9/(2*lamb_temp^2)
    } else if(whiten){
      beta_est[,k] = MASS::ginv(t(W))%*%beta_temp / sqrt(v[k])
    } else {
      beta_est[,k] = beta_temp
    }
    
  }
  
  if((whiten&linear)){
    return_list = list(beta_est=beta_est, pi_est=pi.power)
  } else {
    return_list = list(beta_est=beta_est)
  }
  
  return(return_list)
}



powerMethod_mat = function(X, y, K, GaussianInput=TRUE, linear=FALSE, sparse=FALSE){
  # X is observed predictors
  # y is observed responses
  # K is number of clusters
  
  p = dim(X)[2]
  beta_est = matrix(0,p,K)
  
  T_moment = moment2nd(y,X,linear=linear)
  
  if(sparse){
    SPCtune = SPC.cv(T_moment, sumabsvs=seq(1.1,sqrt(ncol(T_moment)),len=10), trace=FALSE)
    beta_est = SPC(T_moment, sumabsv=SPCtune$bestsumabsv, K=K, trace=FALSE)$v
  } else {
    beta_est = eigen(T_moment)$vectors[,1:K]
  }
  
  return_list = list(beta_est=beta_est, T_moment=T_moment)
  return(return_list)
}


powerMethod_mirror = function(X, y, K, GaussianInput=TRUE, sparse=FALSE){
  # X is observed predictors
  # y is observed responses
  # K is number of clusters
  
  n = dim(X)[1]
  p = dim(X)[2]
  beta_est = matrix(0,p,K)
  
  n1 = floor(n/2)
  d_mirr = rep(0,p)
  for (i in 1:n) {
    d_mirr = d_mirr + y[i]*X[i,]
  }
  d_mirr = d_mirr/n
  
  z = y
  for (i in 1:n) {
    z[i] = y[i] * sign(sum(d_mirr*X[i,]))
  }
  
  
  T_moment = moment2nd(z,X,linear=FALSE)
  
  
  eigen = eigen(T_moment)
  eigen_vals = eigen$values
  eigen_ind = sort(abs(eigen_vals),decreasing=TRUE,index.return=TRUE)$ix
  
  if(sparse){
    SPCtune = SPC.cv(T_moment, sumabsvs=seq(1.2,sqrt(ncol(T_moment)),len=10), trace=FALSE)
    beta_est = SPC(T_moment, sumabsv=SPCtune$bestsumabsv, K=K, trace=FALSE)$v
  } else {
    beta_est = eigen(T_moment)$vectors[,eigen_ind[1:2]]
  }
  
  
  return_list = list(beta_est=beta_est,T_moment=T_moment)
  return(return_list)
}









