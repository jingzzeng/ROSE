# AR function
AR <- function(rho, p){
  m <- matrix(0, p, p)
  for (i in 1:p){
    for (j in 1:p){
      m[i,j] <- rho**(abs(i-j))
    }
  }
  return(m)
}

# CS function
CS <- function(rho, p){
  m <- matrix(rho,p,p)
  diag(m) <- 1
  return(m)
}

# If some whole column is significantly small then we cut the columns to zero
cut_mat <- function(Beta, thrd, rank){
  l <- length(Beta)
  for (i in 1:l){
    if(is.null(Beta[[i]])) next
    mat <- as.matrix(Beta[[i]])
    nobs <- nrow(mat)
    nvars <- ncol(mat)
    r <- rank[i]
    if(r == 0){
      Beta[[i]] <- matrix(0, nobs, nvars)
    }else{
      vec <- as.vector(mat)
      vec[abs(vec) < thrd] <- 0
      Beta[[i]] <- matrix(vec, nobs, nvars)
    }
  }
  return(Beta)
}

# rank function
rank_func <- function(B, thrd){
  d <- svd(B)$d
  r <- sum(d >= thrd)
  return(r)
}

# subspace function for matrices with the same column dimension
subspace <- function(A,B){
  Pa <- qr.Q(qr(A))
  Pa <- Pa %*% t(Pa)
  Pb <- qr.Q(qr(B))
  Pb <- Pb %*% t(Pb)
  d <- dim(A)[2]
  return(norm(Pa-Pb, type="F")/sqrt(2*d))
}

eval_dc <- function(Beta, x, y, yclass = NULL, H = H, categorical = FALSE, original_y = FALSE){
  if(is.null(yclass)){
    if(categorical == FALSE){
      ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/H), na.rm=TRUE))
      yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
    }
    else if(categorical == TRUE){
      yclass <- y
    }
  }
  ydummy <- sapply(unique(yclass)[2:length(unique(yclass))], function(i){
    as.numeric(yclass == i)
  })
  p <- dim(x)[2]
  
  xbar <- colMeans(x)
  xc <- x - t(xbar)[rep(1, NROW(x)),]
  gamma <- apply(xc, 1, function(t){sqrt(mean(t^2))})  # gamma hat.
  weight <- 1/gamma^2
  mu <- c((t(weight) %*% x)/sum(weight))
  w <- (x - t(mu)[rep(1,NROW(x)),])/(gamma %*% t(rep(1,p)))
  
  if(!is.list(Beta)){Beta <- list(Beta)}
  l <- length(Beta)
  result <- sapply(seq_len(l), function(i){
    if(is.null(Beta[[i]])){
      NA
    }else{
      mat <- as.matrix(Beta[[i]])
      if(original_y){
        1 - dcor(w %*% mat, y)
      }else{
        1 - dcor(w %*% mat, ydummy)
      }
    }
  })
  return(result)
}

cut_func <- function(x, ub){
  sign(x) * min(abs(x), ub)
}

# Prediction function for regression problem in real data analysis.
err_mse <- function(x_train, y_train, x_val, y_val){
  x_train <- data.frame(x_train)
  x_val <- data.frame(x_val)
  colnames(x_train) <- paste0('X', seq_len(ncol(x_train)))
  colnames(x_val) <- paste0('X', seq_len(ncol(x_val)))
  data <- cbind(y = y_train, x_train)

  # --------- Huber regression ----------- #
  huber_err <- tryCatch({
    model <- rlm(y~., data = data)
    prediction <- predict(model, newdata = x_val)
    mean((prediction - y_val)^2)
  },error = function(e){
    NA
  })

  # --------- Kernel regression ----------- #
  # ker_err <- tryCatch({
  #   model <- np::npreg(as.formula(paste(c(names(data)[1], paste(names(data)[-1], collapse = '+')), collapse = '~')), data = data, regtype = 'lc', ckernel='gaussian')
  #   prediction <- predict(model, newdata = x_val)
  #   mean((prediction - y_val)^2)
  # },error = function(e){
  #   NA
  # })
  
  # # --------- SVM for regression ----------- #
  # svm_err <- tryCatch({
  #   model <- svm(y~., data = data, scale = FALSE)
  #   prediction <- predict(model, newdata = x_val)
  #   mean((prediction - y_val)^2)
  # },
  # error = function(e){
  #   NA
  # })

  # --------- Random forest ----------- #
  rf_err <- tryCatch({
    model <- randomForest(y~., data = data)
    prediction <- predict(model, newdata = x_val)
    mean((prediction - y_val)^2)
  },
  error = function(e){
    NA
  })

  # --------- Xgboost ----------- #
  xgb_err <- tryCatch({
    xgb_train <- xgb.DMatrix(data = as.matrix(x_train), label = y_train)
    model <- xgb.train(data = xgb_train, nrounds = 3000, objective="reg:squarederror", verbose = 0)
    prediction <- predict(model, newdata = as.matrix(x_val))
    mean((prediction - y_val)^2)
  },
  error = function(e){
    NA
  })
  
  # c(huber_err, ker_err, rf_err, xgb_err)
  c(huber_err, rf_err, xgb_err)
}

dc_score <- function(x, y){
  ybreaks <- as.numeric(quantile(y, probs=seq(0,1, by=1/5), na.rm=TRUE))
  yclass <- cut(y, breaks = ybreaks, include.lowest = TRUE, labels = FALSE)
  ydummy <- sapply(unique(yclass)[2:length(unique(yclass))], function(i){
    as.numeric(yclass == i)
  })
  dcor(ydummy, x)
}


# ## --------------------------------------------------------- ##
# ## Alternative cv criteria: from Tan et al., 2018.
# ## --------------------------------------------------------- ##
# predicty <- function(x, xtrain_new, ytrain){
#   res <- xtrain_new - t(x)[rep(1,NROW(xtrain_new)),]
#   weight <- exp(-0.5*apply((res)^2,1,sum))
#   weight <- weight/sum(weight)
#   sum(weight*ytrain)
# }
# 
# eval_error <- function(Beta, xtrain, ytrain, xtest, ytest){
#   if(!is.list(Beta)){Beta <- list(Beta)}
#   l <- length(Beta)
# 
#   # lb <- quantile(ytrain, 0.1)[[1]]
#   # ub <- quantile(ytrain, 0.9)[[1]]
#   # ytrain <- sapply(ytrain, cut_func, lb = lb, ub = ub)
#   ub <- quantile(abs(ytrain), 0.8)[[1]]
#   ytrain <- sapply(ytrain, cut_func, ub = ub)
# 
#   # lb <- quantile(ytest, 0.1)[[1]]
#   # ub <- quantile(ytest, 0.9)[[1]]
#   # ytest <- sapply(ytest, cut_func, lb = lb, ub = ub)
#   ub <- quantile(abs(ytest), 0.8)[[1]]
#   ytest <- sapply(ytest, cut_func, ub = ub)
# 
#   result <- sapply(seq_len(l), function(i){
#     if(is.null(Beta[[i]])){
#       NA
#     }else{
#       mat <- as.matrix(Beta[[i]])
#       xtrain_new <- xtrain %*% mat
#       xtest_new <- xtest %*% mat
#       # yhat <- apply(X = xtest_new, MARGIN = 1, FUN = predicty, xtrain_new, ytrain)
#       yhat <- rep(0, NROW(xtest_new))
#       for(i in 1:NROW(xtest_new)){
#         yhat[i] <- predicty(xtest_new[i,], xtrain_new, ytrain)
#       }
#       error <- mean((ytest - yhat)^2)
#       error
#     }
#   })
# 
#   result
# }
