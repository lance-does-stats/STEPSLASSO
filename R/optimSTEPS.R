#' Optimize the likelihood function values for high-dimensional sparse regression for secondary trait analysis under extreme phenotype sampling data.
#'
#' This function optimizes the likelihood function for high-dimensional extreme phenotype sampling data.
#' @param data.mat Data frame containing Y, Z, X, c1, c2. Required.
#' @param c1 Right censored point. Required.
#' @param c2 Left censored point. Required.
#' @import MASS magrittr stats
#' @export

optimSTEPS <- function(data.mat){


  c1 <- data.mat$c1
  c2 <- data.mat$c2
  df.optim <- data.mat$data

  Y <- df.optim[,1]
  X <- df.optim[,-c(1,2)] %>% as.matrix()
  Z <- df.optim[,2]

  n <- nrow(X)
  K <- ncol(X)

  if(dim(X)[2]==0){
    n <- length(Y)
    K <- 2 #done to pass an empty matrix as covariates; K=1 vector messes with
    X <- matrix(0, nrow=n , ncol= K)

    tag.noX = TRUE
  } else tag.noX = FALSE


  if(dim(X)[2]==1){
    n <- length(X)
    K <- 2
    X <- cbind(X,0)

    tag.oneX = TRUE
  } else tag.oneX = FALSE



  if(n>K){
    alpha_0 = as.vector(summary(lm(Z~X))$coefficients[,1])
    beta_0 = as.vector(summary(lm(Y~Z+X))$coefficients[-2,1])

    gam_0 = as.vector(summary(lm(Y~Z+X))$coefficients[2,1])
    sd.z_0 = as.vector(summary(lm(Z~X))$sigma)
    sd.y_0 <- data.mat$sd.y

  }else{
    cv.fit.Z = cv.glmnet(X,Z, alpha=0, standardize = FALSE)
    bestlam.Z     = cv.fit.Z$lambda.min
    temp.Z     = glmnet(X,Z,lambda=bestlam.Z, standardize = FALSE, alpha=0)
    alpha_0 <- c(as.vector(temp.Z$a0),as.vector(temp.Z$beta))

    cv.fit.Y = cv.glmnet(X,Y,standardize = FALSE, alpha=0)
    bestlam.Y     = cv.fit.Y$lambda.min
    temp.Y     = glmnet(X,Y,lambda=bestlam.Y, standardize = FALSE, alpha=0)
    beta_0 <- c(as.vector(temp.Y$a0),as.vector(temp.Y$beta))

    gam_0 = as.vector(summary(lm(Y~Z))$coefficients[2,1])
    sd.z_0 = (K+n+1)/(n^2 +n)*norm(Z, type="2")^2 - norm(t(X)%*%Z,type="2")^2/(n^2 +n)
    sd.y_0 <- data.mat$sd.y
  }



  if(tag.noX == TRUE){
    alpha_0 <- c(alpha_0, rep(0,K))
    beta_0 <- c(beta_0, rep(0,K))
  }

  if(tag.oneX == TRUE){
    alpha_0 <- c(alpha_0, rep(0,K-1))
    beta_0 <- c(beta_0, rep(0,K-1))
  }

df.optim2 <- cbind(Y,Z,X) %>% as.matrix()

par.in <- c(alpha_0, beta_0, gam_0, sd.z_0, sd.y_0)


fn = function(par){
  para.ll = list(A0= par[1], A = par[1+1:K],
                 B0=par[2+K], B = par[2+K+1:K],
                 gam = par[3+2*K],
                 sd.z = par[4+2*K],
                 sd.y = par[5+2*K])
  res.ll = likelihoodSTEPS(df.optim2, c1=c1, c2=c2, para.ll)
  return(res.ll$ll)
}

gr = function(par){
  para.ll = list(A0= par[1], A = par[1+1:K],
                 B0=par[2+K], B = par[2+K+1:K],
                 gam = par[3+2*K],
                 sd.z = par[4+2*K],
                 sd.y = par[5+2*K])
  res.ll = likelihoodSTEPS(df.optim2, c1=c1, c2=c2, para.ll)
  return(res.ll$d.ll)
}
optimized=1

res.opt=tryCatch(optim(par.in, fn, gr, method = "L-BFGS-B",
                       lower=c(-10000, rep(-10000,K), -10000, rep(-10000,K), -1000, 0.0001, 0.0001),
                       upper=c(10000, rep(10000,K), 10000, rep(10000,K), 1000, 1000, 1000)),
                 error=function(e) e)

if(is(res.opt,"error")){
  res.opt=tryCatch(optim(par.in,fn,gr,method = "CG"),
                   error=function(e) e)
  method2error=1
}

if(is(res.opt,"error") && method2error==1){
  res.opt=tryCatch(optim(par.in,fn,gr,method = "BFGS"),
                   error=function(e) e)
  method3error=1
}


if(is(res.opt,"error") && method3error==1){
  res.opt=list(par=par.in)
  optimized=0
}


names(res.opt$par) = c("a0", paste("a",1:K,sep=""),
                       "b0", paste("b",1:K,sep=""),
                       "gamma","sdz","sdy")


alpha= res.opt$par[1+1:K]
beta=res.opt$par[2+K+1:K]

if(tag.noX == TRUE){
  alpha <- NA
  beta <- NA
}

if(tag.oneX == TRUE){
  alpha <- alpha[1]
  beta <- beta[1]
}


return(list(alpha0= res.opt$par[1],
            alpha= alpha,
            beta0= res.opt$par[2+K],
            beta= beta,
            gamma=res.opt$par[3+2*K],
            sdz=res.opt$par[4+2*K],
            sdy=res.opt$par[5+2*K],
            optimized=optimized))

}







