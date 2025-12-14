library(MASS)
library(rTensor)
library(combinat)

####################################
##### Subspace Distance metric #####
####################################
subspace = function(A,B){
  # calculate the distance between two subspaces span(A) and span(B)
  # A and B are both p-by-u matrix
  # returns a distance metric that is between 0 and 1
  Pa = qr.Q(qr(A))
  Pa = Pa %*% t(Pa)
  Pb = qr.Q(qr(B))
  Pb = Pb %*% t(Pb)
  u = dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*u))
}


##############################################
##### Median bootstrap to standard error #####
##############################################
med_se = function(x,B){
  B_median = rep(0,B)
  n = length(x)
  for (i in 1:B) {
    id = sample(1:n,n,replace=T)
    B_median[i] = median(x[id])
  }
  return(sd(B_median)/sqrt(B))
}


##########################
##### Outer product ######
##########################
outer3d = function(x,y,z){
  return(outer(x,outer(y,z)))
}



#######################
##### Match norm ######
#######################
norm_match = function(beta,beta_target){
  
  K = dim(beta)[2]
  ind = rep(0,K)
  for (i in 1:K) {
    spc_err = rep(0,K)
    beta_temp = beta[,i]
    for (j in 1:K) {
      spc_err[j] = subspace(matrix(beta_temp,ncol=1), matrix(beta_target[,j],ncol=1))
    }
    ind[i] = which.min(spc_err)
  }
  
  beta_target_new = beta_target[,ind]
  
  norm_target = apply(beta_target_new, 2, norm, "2")
  norm_now = apply(beta, 2, norm, "2")
  
  scalar = norm_target/norm_now
  beta_new = t(t(beta)*scalar)
  
  return_list = list(beta_new=beta_new)
  return(return_list)
}



#######################################################
##### Empirical estimation of second order moment #####
#######################################################
moment2nd = function(y, X, GaussianInput=TRUE, linear=FALSE){
  # y is a vector of sample responses
  # X is a matrix of sample covariates
  n = dim(X)[1]
  p = dim(X)[2]
  Ip = diag(1,p)
  
  if(GaussianInput){
    M = 0
    for (i in 1:n) {
      if(linear==TRUE){
        y_temp = y[i]^2
      }
      else{
        y_temp = y[i]
      }
      
      M = M + y_temp*(outer(X[i,],X[i,])-Ip)
      #M = M + y_temp*(outer(X[i,],X[i,]))
      #M = M + outer(X[i,],X[i,])
      #M = M + outer(X[i,],X[i,]) %*% outer(X[i,],X[i,])
    }
    
    M = M/n
  }
  
  return(M)
}



######################################################
##### Empirical estimation of third order moment #####
######################################################
moment3rd = function(y, X, GaussianInput=TRUE, linear=FALSE){
  # y is a vector of sample responses
  # X is a matrix of sample covariates
  n = dim(X)[1]
  p = dim(X)[2]
  Ip = diag(1,p)
  
  if(GaussianInput){
    M = 0
    for (i in 1:n) {
      if(linear==TRUE){
        y_temp = y[i]^3
      }
      else{
        y_temp = y[i]
      }
      
      temp1 = y_temp*outer3d(X[i,], X[i,], X[i,])
      temp2 = temp3 = temp4 = 0
      for (j in 1:p) {
        temp2 = temp2 + y_temp*outer3d(Ip[j,],X[i,],Ip[j,])
        temp3 = temp3 + y_temp*outer3d(Ip[j,],Ip[j,],X[i,])
        temp4 = temp4 + y_temp*outer3d(X[i,],Ip[j,],Ip[j,])
      }
      temp = temp1 - (temp2+temp3+temp4)
      M = M + temp
    }
    
    M = M/n
  }
  
  return(M)
}



############################
##### Tensor whitening #####
############################
whiten = function(T,r,V=NULL){
  p = dim(T)[1]
  
  if(is.null(V)==FALSE){

    eig = eigen(V)
    U = eig$vectors[,1:r]
    v = eig$values[1:r]
    
    degen = sum(v<0)
    
    if(degen>0) print("degenerated!")
    
  } else {
    degen = 1
  }
  
  mi=1
  while(degen>0){
    
    theta = rnorm(p)
    
    V = 0
    for (i in 1:p) {
      V = V + theta[i]*T[,,i]
    }
    
    V = rTensor::ttm(rTensor::as.tensor(T), matrix(theta,nrow=1), m=3)@data[,,1] 
    eig = eigen(V)
    U = eig$vectors[,1:r]
    v = eig$values[1:r]
    
    degen = sum(v<0)
    mi=mi+1
    if(mi>100){break}
  }
  
  if(degen>0) stop("degenerated error")
  
  
  W = U%*%diag(1/sqrt(v))
  T_whitened = rTensor::ttl(rTensor::as.tensor(T), list(t(W),t(W),t(W)), ms=c(1:3))@data
  
  return_list = list(W=W, T_whitened=T_whitened, V=V, v=v)
  return(return_list)
}
 


####################################
##### SVD-based initialization #####
####################################
svd_init = function(T){
  r = dim(T)[1]
  theta = rnorm(r)
  
  V = rTensor::ttm(rTensor::as.tensor(T), matrix(theta,nrow=1), m=3)@data[,,1]
  
  u = eigen(V)$vectors[,1]
  
  return_list = list(u=u)
  return(return_list)
}



#############################################
##### Rank-1 tensor power decomposition #####
#############################################

# powerDecomp_rank1 = function(T,N,L,nu){
#   p = dim(T)[1]
#   Theta = matrix(0,p,L)
#   for (i in 1:L) {
#     theta_now = svd_init(T)$u
#     for (t in 1:N) {
#       mat_list = list(diag(1,p),matrix(theta_now,nrow=1),matrix(theta_now,nrow=1))
#       theta_now = rTensor::ttl(rTensor::as.tensor(T), mat_list, ms=1:3)@data
#       theta_now = as.vector(theta_now)
#       theta_now = theta_now / norm(theta_now,type="2")
#     }
#     Theta[,i] = theta_now
#   }
#   
#   iter = 1
#   while((ncol(Theta)>0)&(iter<100)){
#     s = ncol(Theta)
#     
#     lamb_list = rep(0,s)
#     for (i in 1:s) {
#       mat_list = list(matrix(Theta[,i],nrow=1),matrix(Theta[,i],nrow=1),matrix(Theta[,i],nrow=1))
#       lamb_list[i] = abs(rTensor::ttl(rTensor::as.tensor(T), mat_list, ms=1:3)@data[1])
#     }
#     ind_max = which.max(lamb_list)
#     
#     theta_now = Theta[,ind_max]
#     for (t in 1:N) {
#       mat_list = list(diag(1,p),matrix(theta_now,nrow=1),matrix(theta_now,nrow=1))
#       theta_now = rTensor::ttl(rTensor::as.tensor(T), mat_list, ms=1:3)@data
#       theta_now = as.vector(theta_now)
#       theta_now = theta_now / norm(theta_now,type="2")
#     }
#     
#     inner_list = rep(0,s)
#     for (i in 1:s) {
#       inner_list[i] = abs(sum(Theta[,i]*theta_now))
#     }
#     
#     Theta = Theta[,!(inner_list>nu/2),drop=FALSE]
#     
#     iter = iter+1 
#     
#     #print(s)
#     #print(iter)
#   }
#   
#   
#   theta = theta_now
#   mat_list = list(matrix(theta,nrow=1),matrix(theta,nrow=1),matrix(theta,nrow=1))
#   lamb = rTensor::ttl(rTensor::as.tensor(T), mat_list, ms=1:3)@data[1]
#   
#   return_list = list(lamb=lamb, theta=theta)
#   return(return_list)
# }




##############################################################
##### Rank-1 tensor power decomposition from Zeng, et al #####
##############################################################

powerDecomp_rank1 = function(T,N,L){
  p = dim(T)[1]
  Theta = matrix(0,p,L)
  for (i in 1:L) {
    theta_now = svd_init(T)$u
    for (t in 1:N) {
      mat_list = list(diag(1,p),matrix(theta_now,nrow=1),matrix(theta_now,nrow=1))
      theta_now = rTensor::ttl(rTensor::as.tensor(T), mat_list, ms=1:3)@data
      theta_now = as.vector(theta_now)
      theta_now = theta_now / norm(theta_now,type="2")
    }
    Theta[,i] = theta_now
  }

  lamb_list = rep(0,L)
  for (i in 1:L) {
    mat_list = list(matrix(Theta[,i],nrow=1),matrix(Theta[,i],nrow=1),matrix(Theta[,i],nrow=1))
    lamb_list[i] = abs(rTensor::ttl(rTensor::as.tensor(T), mat_list, ms=1:3)@data[1])
  }
  ind_max = which.max(lamb_list)
  
  theta_now = Theta[,ind_max]
  for (t in 1:N) {
    mat_list = list(diag(1,p),matrix(theta_now,nrow=1),matrix(theta_now,nrow=1))
    theta_now = rTensor::ttl(rTensor::as.tensor(T), mat_list, ms=1:3)@data
    theta_now = as.vector(theta_now)
    theta_now = theta_now / norm(theta_now,type="2")
  }


  theta = theta_now
  mat_list = list(matrix(theta,nrow=1),matrix(theta,nrow=1),matrix(theta,nrow=1))
  lamb = rTensor::ttl(rTensor::as.tensor(T), mat_list, ms=1:3)@data[1]

  return_list = list(lamb=lamb, theta=theta)
  return(return_list)
}










