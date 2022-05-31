#' Perform Wald test for ST-EPS data
#'
#' @param data.mat List containing Y, Z, X, c1, c2. Required.
#' @param position indicates position in coefficient matrix of variable to test. Required.
#' @param null_value the value to test the null hypothesis against (default=0)
#' @import stats
#' @export
#' @example
#' out_score <- lrtTestSTEPS(data.mat=data.mat, position=1, null_value=0)


lrtTestSTEPS <- function(data.mat, position=1){

  df <- data.mat$data
  X <- df[,-c(1,2)] %>% as.matrix()
  Z <- df[,2]


  if(dim(X)[2]==1){
    full.model <- stats::glm(Z~X)
    null.model <- stats::glm(Z~1)
  }else{
    full.model <- stats::glm(Z~X)
    null.model <- stats::glm(Z~X[,-position])
  }

  # LR test
  p.value.LRT <- stats::anova(null.model,full.model,test="LRT")[2,5]
  Z.lrt <- qnorm(p.value.LRT)

  list('LRT P-value'=p.value.LRT, 'LRT Z-value'=Z.lrt)

}
