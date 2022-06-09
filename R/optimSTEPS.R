#' Optimize the likelihood function values for high-dimensional sparse regression for secondary trait analysis under extreme phenotype sampling data.
#'
#' This function optimizes the likelihood function for high-dimensional extreme phenotype sampling data.
#' @param data.mat Data frame containing Y, Z, X, c1, c2. Required.
#' @param c1 Right censored point. Required.
#' @param c2 Left censored point. Required.
#' @import MASS magrittr stats
#' @export

optimSTEPS <- function(data.mat, method="L-BFGS-B"){


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



alpha_0 = as.vector(summary(lm(Z~X))$coefficients[,1])
beta_0 = as.vector(summary(lm(Y~Z+X))$coefficients[-2,1])
gam_0 = as.vector(summary(lm(Y~Z+X))$coefficients[2,1])
sd.z_0 = as.vector(summary(lm(Z~X))$sigma)
sd.y_0 = as.vector(summary(lm(Y~X+Z))$sigma)

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

if(method=="CG" || method == "BFGS"){
  res.opt=optim(par.in,fn,gr,method = method)
} else{
  res.opt=optim(par.in, fn, gr, method = method,
                lower=c(0, rep(0,pX), 0, rep(0,pX), 0, 0.0001, 0.0001),
                upper=c(10000, rep(10000,pX), 10000, rep(10000,pX), 1000, 1000, 1000))
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
            sdy=res.opt$par[5+2*K]))

}







