# Generate ST-EPS data
# param: n - number of individuals in full sample to generate
# param: p - number of parameters to generate (MVN, iid)
# param: q.a - number of nonzero covariates for secondary trait
# paragm: q.b - number of nonzero covariates for primary trait
# param: A - px1 matrix of true secondary trait parameter values
# param: B - px1 matrix of true primary trait paramter values
# param: gam - true correlation between primary and secondary traits (default=1)
# param: sd.y - true standard deviation of primary trait (default=1)
# param: sd.z - true standard deviation of secondary trait (default=1)
# param: eps - upper and lower quantile to sampled from full sample to phenotype (default=0.5; full sample)

# required libraries
#library(MASS)
#library(dplyr)


generateSTEPS <- function(n, p, q.a, q.b, A_int=0, B_int=0, A_coef, B_coef, overlap=0, gam=1, sd.y=1, sd.z=1, eps=0.5){
  
  A <- c(A_int,rep(A_coef,q.a),
          rep(0,p-q.a))
  
  B <- c(B_int,rep(0,q.a-overlap),
          rep(B_coef,q.b),
          rep(0,p-q.b-q.a+overlap))
  
  if(length(A)!=p+1) stop('Wrong dimensions of secondary trait coefficients.')
  if(length(B)!=p+1) stop('Wrong dimensions of primary trait coefficients.')
  if(sd.y<=0 || sd.z<=0) stop('Standard deviation must be positive, non-zero.')
  
  X <- MASS::mvrnorm(n = n, mu = rep(0,p), Sigma = diag(p)) %>% base::as.matrix()
  
  Z <- A[1] + X%*%A[-1] + rnorm(n, 0, sd.z) %>% base::as.matrix()
  
  Y <- B[1] + X%*%B[-1] + gam*Z + rnorm(n, 0, sd.y) %>% base::as.matrix()
  
  temp <- data.frame(Y,Z,X)
  
  c1 <- quantile(Y,eps)[[1]]
  c2 <- quantile(Y,1-eps)[[1]]
  
  ok <-  which(temp$Y>c2|temp$Y<c1)
  #('c1=', c1, '\n')
  #cat('c2=', c2, '\n')
  data <- temp[ok,] %>% base::as.matrix()
  
  list("A0"=A[1], "A"=A[-1], "B0"=B[1], "B"=B[-1], "c1"=c1, "c2"=c2, "gam"=gam, "sd.y"=sd.y, "sd.z"=sd.z, "data"=data)
  
}

#data.mat <- generateSTEPS(n=1000, p=10, q.a=3, q.b=3, A_coef = 0.4, B_coef = 0.4, gam=0.8, sd.y=1, sd.z=1, eps=0.5)

