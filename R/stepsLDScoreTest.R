#' Perform Score test for ST-EPS data
#'
#' @param data.mat List containing Y, Z, X, c1, c2. Required.
#' @param position indicates position in coefficient matrix of variable to test. Required.
#' @param null_value the value to test the null hypothesis against (default=0)
#' @import
#' @export
#' @example
#' out_score <- scoreTestSTEPS(data.mat=data.mat, position=1, null_value=0)


scoreTestSTEPS <- function(data.mat, position=1, null_value=0){

  c1 <- data.mat$c1
  c2 <- data.mat$c2
  df <- data.mat$data


  Y <- df[,1]
  X <- df[,-c(1,2)]
  Z <- df[,2]

  estimates <- optimSTEPS(data.mat)
  theta <- sqrt(estimates$gamma^2 * estimates$sdz^2 + estimates$sdy^2) %>% as.numeric()
  A0 <- estimates$alpha0
  A <- estimates$alpha
  As <- A; As[position] <- null_value
  B0 <- estimates$beta0
  B <- estimates$beta
  gam <- estimates$gamma
  sdz <- estimates$sdz
  sdy <- estimates$sdy

  if(length(A)==1){
    n <-  length(X)
    K <-  1
    XA <- A0 + X*As
    XB <- B0 + X*B
  }else{
    n <-  nrow(X)
    K <-  ncol(X)
    XA <- A0 + X%*%As
    XB <- B0 + X%*%B
  }

  t1 <- c1 - (XB + gam*XA)# lower
  t2 <- c2 - (XB + gam*XA)# upper

  f1 <- dnorm(t1, sd=theta)
  f2 <- dnorm(t2, sd=theta)
  F1 <- pnorm(t1, sd=theta)
  F2 <- pnorm(t2, sd=theta)
  P.S1 <- F1 + 1 - F2

  mv=(f1-f2)/P.S1
  vv=(t1*f1-t2*f2)/P.S1 + mv^2

  # Score function
  score.alpha <- t(X)%*%(((Z-XA)/sdz^2) + (f1 - f2)/P.S1*gam/theta)

  # Fisher's information matrix
  J = matrix(0,nrow=K,ncol=K)
  J[1:K, 1:K]<-t(X) %*%diag(as.numeric(-1/sdz^2 + (gam^2 / theta^2)*(vv)))%*%X
  fish.info <- -J

  # Score test
  score.value <- (score.alpha[position])^2 / fish.info[position,position]
  Z.s <- score.alpha[position] / sqrt(fish.info[position,position])
  p.value.score <- 1 - pchisq(score.value, df=1)


  list('Estimate'=A[position],'Score P-value'=p.value.score, 'Score Z-value'=Z.s)

}

