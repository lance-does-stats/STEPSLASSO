#' Calculate initial values for sdz and gamma using GLMNET for high-dimensional sparse regression for secondary trait analysis under extreme phenotype sampling data.
#'
#' This function calculates initial values for sdz and gamma for high-dimensional extreme phenotype sampling data.
#' @param data.mat Data frame containing Z, X. Required.
#' @import magrittr glmnet
#' @export


stepsGLMNET <- function(data.mat){

  ### Step 1a) get covariates to include

  # split data into training and testing data to use CV for best lambda
  train <- sample(1:nrow(data.mat$data[,-c(1,2)]), nrow(data.mat$data[,-c(1,2)])/2)
  test <- (-train)
  z.test <- data.mat$data[test,2]

  #Use CV to choose best lambda based on MSE
  cvfit <- glmnet::cv.glmnet(x=data.mat$data[train,-c(1,2)],y=data.mat$data[train,2], alpha=1)
  bestlam <- cvfit$lambda.min

  grid <- 10^ seq (10 , -2 , length = 100)
  out <- glmnet::glmnet(x=data.mat$data[,-c(1,2)], y=data.mat$data[,2], alpha=1, lambda=grid)


  alpha.k0 <- stats::predict(out, type="coefficients", s=bestlam) %>% as.vector()
  include <- which(alpha.k0!=0)[-1] + 1


  data.mat$data <- data.mat$data[,c(1,2,include)]


  ### Step 1b) get estimates for sdz, sdy, gamma
  if(length(include)>=dim(data.mat$data)[1]){
    X.temp <- data.mat$data[,-c(1,2)]
    Y.temp <- data.mat$data[,1]
    Z.temp <- data.mat$data[,2]
    K <- ncol(X.temp)
    n <- nrow(X.temp)
    Sigma.hat <- cov(X.temp)
    sdz <- (K+n+1)/(n^2 +n)*norm(Z.temp, type="2")^2 - norm(Sigma.hat%*%t(X.temp)%*%Z.temp,type="2")^2/(n^2 +n)
    #sdz <- (K+n+1)/(n^2 +n)*norm(Z.temp, type="2")^2 - norm(t(X.temp)%*%Z.temp,type="2")^2/(n^2 +n)
    gam <- lm(Y.temp~Z.temp)$coef[2] %>% as.numeric()
  }else{
    estimates <- optimSTEPS(data.mat)
    gam <- estimates$gamma
    sdz <- estimates$sdz
  }


  list(int.sdz=sdz, int.gamma=gam)

}
