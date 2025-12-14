AR <- function(rho, p){
  Sigma <- diag(p)
  for (i in 1:p){
    for (j in 1:p){
      Sigma[i,j] <- rho^abs(i-j)
    }
  }
  return(Sigma)
}

ERM <- function(p){
  omega1 <- matrix(rbinom(p*p,1,0.1),p,p)
  omega2 <- matrix(runif(p*p,0.5,1),p,p)
  omega3 <- sign(matrix(rnorm(p*p),p,p))
  omega4 <- omega1*omega2*omega3
  omega5 <- (omega4+t(omega4))/2
  omega6 <- (omega5)+(max(-eigen(omega5)$values[p],0)+0.05)*diag(p)
  w <- 1/as.vector(sqrt(diag(omega6)))
  wwt <- w%*%t(w)
  omega7 <- omega6*wwt
  return(omega7)
}

predict.test <- function(x, y, beta, w, sigma, K){
  tmp <- y %*% t(rep(1,K)) - x %*% beta
  D <- diag(w)
  post <- dnorm(tmp, mean = 0, sd = sigma)
  num <- post %*% D
  denom <- apply(num, 1, sum)
  eta <- num/(denom %*% t(rep(1,K)))
  if(any(denom == 0)) eta[which(denom == 0),] <- w
  prediction <- apply(eta, 1, which.max)
  prediction
}

update.sigma <- function(x, y, beta, eta){
  n <- NROW(x)
  K <- NCOL(beta)
  tmp <- (y %*% t(rep(1,K)) - x %*% beta)^2
  sigma <- sum(tmp * eta)/n
  sigma
}