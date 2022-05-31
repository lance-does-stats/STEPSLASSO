#' Calculate output for stepsLASSO for high-dimensional sparse regression for secondary trait analysis under extreme phenotype sampling data.
#'
#' This function calculates final refit estimates and pvalues for high-dimensional extreme phenotype sampling data.
#' @param data.mat Data frame containing Y,Z,c1,c2, and selected covariates X from stepsLassoSolver. Required.
#' @import
#' @export


stepsLD <- function(data.mat){
  Y <- data.mat$data[,1]
  Z <- data.mat$data[,2]
  X <- data.mat$data[,-c(1,2)]

  estimates <- optimSTEPS(data.mat)

  #number of selected covariates - k
  k <- length(estimates$alpha)

  # Initiate vectors
  pval.s <- c()
  pval.w <- c()
  pval.lrt <- c()

  # Calculate 3 types of p-value for every selected covariate
  if(k==1){
    pval.s <- scoreTestSTEPS(data.mat, position=1)$`Score P-value`
    pval.w <- waldTestSTEPS(data.mat, position=1)$`Wald P-value`
    pval.lrt <- lrtTestSTEPS(data.mat, position=1)$`LRT P-value`
  }else{
    for(i in 1:k){
      pval.s[i] <- scoreTestSTEPS(data.mat, position=i)$`Score P-value`
      pval.w[i] <- waldTestSTEPS(data.mat, position=i)$`Wald P-value`
      pval.lrt[i] <- lrtTestSTEPS(data.mat, position=i)$`LRT P-value`
    }
  }

  # list output
  list(alpha.table=cbind(estimates$alpha, pval.s, pval.w, pval.lrt),
       alpha.hat = estimates$alpha,
       sdz.hat=estimates$sdz,
       gamma.hat=estimates$gamma,
       sdy.hat=estimates$sdy,
       beta.hat=estimates$beta)

}



