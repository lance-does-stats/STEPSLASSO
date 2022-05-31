#' Calculate likelihood function values for high-dimensional sparse regression for secondary trait analysis under extreme phenotype sampling data.
#'
#' This function calculates the likelihood function values for high-dimensional extreme phenotype sampling data.
#' @param data.mat Data frame containing Y, Z, X. Required.
#' @param c1 Right censored point. Required.
#' @param c2 Left censored point. Required.
#' @param para Vector of parameter estimates. Required.
#' @import MASS magrittr
#' @export
#' @example
#' df <- generateSTEPS(n=1000, p=10, q.a=3, q.b=3, A_coef=0.8, B_coef=0.8, gam=0.8, sd.y=1, sd.z=1, eps=0.2)
#' c1 <- df$c1
#' c2 <- df$c2
#' data.mat <- df$data
#' para <- df[-10]
#' likelihoodSTEPS(df2, c1=c1, c2=c2, para = para.in)

likelihoodSTEPS <- function(data.mat, c1, c2, para){

  X <- data.mat[,-c(1,2)]
  Y <- data.mat[,1]
  Z <- data.mat[,2]
  A0 <- para$A0
  A <- para$A
  B0 <- para$B0
  B <- para$B
  gam <- para$gam
  sdy <- para$sd.y
  sdz <- para$sd.z
  theta <- sqrt(gam^2 * sdz^2 + sdy^2)

  n <-  nrow(X)
  K <-  ncol(X)

  XA <- A0 + X%*%A
  XB <- B0 + X%*%B

  t1 <- c1 - (XB + gam*XA) # lower
  t2 <- c2 - (XB + gam*XA)  # upper
  F1 <- pnorm(t1, sd=theta)
  F2 <- pnorm(t2, sd=theta)
  P.S1 <- F1 + 1 - F2

  # Joint distribution of Y & Z
  f <- as.vector(dnorm(Z-XA, sd=sdz)*dnorm(Y-XB-gam*Z, sd=sdy) / P.S1)

  f1 <- dnorm(t1, sd=theta)
  f2 <- dnorm(t2, sd=theta)
  z..XA <- Z-XA
  y..XB..gZ <- Y - XB - gam*Z

  # intermediate step to derivative of f over parameter c(alpha,beta,theta,gamma,sdz,sdy)
  df.int <- f*cbind(z..XA/sdz^2 + (f1-f2)/P.S1*gam, #d.alpha
                    y..XB..gZ/sdy^2 + (f1-f2)/P.S1, #d.beta
                    (t1*f1 - t2*f2) / P.S1 / theta, #d.theta
                    (y..XB..gZ/sdy^2)*Z + ((f1-f2) / P.S1)*X%*%A, #d.gam
                    (z..XA^2 - sdz^2)/sdz^3, #d.sdz
                    (y..XB..gZ^2 - sdy^2)/sdy^3#d.sdy
  )

  a <- cbind(1,X)
  # derivative of f over parameter c(alpha,beta,theta,gamma,sdz,sdy)
  df <- cbind(df.int[,1]*a, #a1, a2, ...
              df.int[,2]*a, #b1, b2, ...
              df.int[,3]*gam*sdz^2 / theta + df.int[,4], #gamma
              df.int[,3]*(gam^2)*sdz / theta + df.int[,5], #sd.z
              df.int[,3]*sdy / theta + df.int[,6] #sd.y
              )

  colnames(df) = c("a0", paste("a",1:K,sep=""),
                   "b0", paste("b",1:K,sep=""),
                   "gamma","sd.z","sd.y")

  ll <- mean(-1*log(f))
  d.ll <- colMeans(-1*df/f)

  return(list(ll=ll, d.ll=d.ll))

}









