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
    pval.s.int <- tryCatch(scoreTestSTEPS(data.mat, position=1), error=function(e) e)
    if(is(pval.s.int,"error")){
      pval.s= -1
    } else{
      pval.s <- pval.s.int$`Score P-value`
    }


    pval.w.int <- tryCatch(waldTestSTEPS(data.mat, position=1), error=function(e) e)
      if(is(pval.w.int,"error")){
        pval.w= -1
      } else{
        pval.w <- pval.w.int$`Wald P-value`
      }

    pval.lrt.int <- tryCatch(lrtTestSTEPS(data.mat, position=1), error=function(e) e)
    if(is(pval.lrt.int,"error")){
      pval.lrt= -1
    } else{
      pval.lrt <- pval.lrt.int$`LRT P-value`
    }


  }else{
    for(i in 1:k){
      pval.s.int <- tryCatch(scoreTestSTEPS(data.mat, position=i), error=function(e) e)
      if(is(pval.s.int,"error")){
        pval.s[i]= -1
      } else{
        pval.s[i] <- pval.s.int$`Score P-value`
      }

      pval.w.int <- tryCatch(waldTestSTEPS(data.mat, position=i), error=function(e) e)
      if(is(pval.w.int,"error")){
        pval.w[i]= -1
      } else{
        pval.w[i] <- pval.w.int$`Wald P-value`
      }

      pval.lrt.int <- tryCatch(lrtTestSTEPS(data.mat, position=i), error=function(e) e)
      if(is(pval.lrt.int,"error")){
        pval.lrt[i]= -1
      } else{
        pval.lrt[i] <- pval.lrt.int$`LRT P-value`
      }

    }
  }

  # list output
  list(alpha.table=cbind(estimates$alpha, pval.s, pval.w, pval.lrt),
       alpha.hat = estimates$alpha,
       sdz.hat=estimates$sdz,
       gamma.hat=estimates$gamma,
       optimized=estimates$optimized)

}



