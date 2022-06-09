#' Test for high-dimensional sparse regression for secondary trait analysis under extreme phenotype sampling data.
#'
#' This function tests the effects of predictors modeled jointly for high-dimensional extreme phenotype sampling data.
#' @param Y Vector of primary trait values. Required.
#' @param c1 Left censored point. Required.
#' @param c2 Right censored point. Required.
#' @param Z Vector of secondary trait values. Required.
#' @param X Matrix of predictors. Required. Required.
#' @param beta.hat Vector of estimated coefficients for eps(Y)~Z. Required.
#' @param sdy.hat Numerical value of estimated variance for eps(Y)~Z. Required.
#' @param maxIter The maximum refining steps when m_w="dzg". Default is 1000.
#' @param verbose Print debugging info or not.
#' @import parallel doParallel foreach
#' @export
#' @examples
#' data.mat0 <- generateSTEPS(n=1000, p=1000, q.a=3, q.b=3, A_coef=0.1, B_coef=0.1, overlap=0, gam=0.12, sd.y=1, sd.z=1, eps=0.2)
#' c1 <- data.mat0$c1
#' c2 <- data.mat0$c2
#' df.0 <- data.mat0$data
#' Y <- df.0[,1]
#' beta.hat <- data.mat0$B
#' sdy.hat <- data.mat0$sd.y
#' X <- df.0[,-c(1,2)]
#' Z <- df.0[,2]
#' stepsLasso1(Y=Y, c1=c1, c2=c2, Z=Z, X=X, beta.hat = beta.hat, sdy.hat = sdy.hat)


stepsLasso <- function(Y, c1, c2, Z, X, beta.hat, sdy.hat, maxIter=1000, verbose=FALSE){
  #### STEP 1 - Estimate Beta.hat and sdy.hat from EPS-LASSO ####

  beta.hat <- beta.hat
  sdy.hat <- sdy.hat



  #### STEP 2 - Get initial estimate for sdz and gamma from GLMNET ####

  data.mat2 <- list(c1=c1, c2=c2, data=cbind(Y,Z,X))
  estimates2 <- stepsGLMNET(data.mat2)

  # initial estimate for sdz
  sdz2 <- estimates2$int.sdz
  gam2 <- estimates2$int.gamma


  #### STEP 3 - Choose best lambda to pass to stepsLassoSolver ####

  nX <-  nrow(X)
  pX <-  ncol(X)


  # Create range of lambda for optimization
  lam0=NULL
  if(is.null(lam0)){
    lamMax=max(abs(t(X)%*%Z))
    if(nX<pX){
      lamMin=0.01*lamMax
    }else{
      lamMin=0.0001*lamMax
    }
    lamR=exp((log(lamMax)-log(lamMin))/99)
    lam0=c(lamMax,lamMax*lamR^(-seq(1:99)))
    lam0=lam0[order(lam0)];
  }

  df3 <- cbind(Y,Z,X) %>% as.matrix()
  #### HDIC ####
  # fn = function(par){
  #   para.ll = list(A0= par[1], A = par[1+1:K],
  #                  B0=par[2+K], B = par[2+K+1:K],
  #                  gam = par[3+2*K],
  #                  sd.z = par[4+2*K],
  #                  sd.y = par[5+2*K])
  #   res.ll = likelihoodSTEPS(df3, c1=c1, c2=c2, para.ll)
  #   return(res.ll$ll)
  # }
  #
  # # Parallel
  # options(warn = -1)
  # cl <- parallel::makeCluster(parallel::detectCores() - 1)
  # doParallel::registerDoParallel(cl)
  # HDIC.out=foreach(i = 1:length(lam0), .combine=c, .export=c('stepsLassoSolver','stepsLeastR','stepsGLMNET','likelihoodSTEPS'), .packages=c('MASS','glmnet')) %dopar%{
  #   HDIC <- c()
  #
  #
  #   iter <- stepsLassoSolver(A=X, Y2=Z, Y1=Y, X1=beta.hat, gamma=gam2, c1=c1, c2=c2,
  #                            lambda = lam0[i], sigma2=sdz2, sigma1 = sdy.hat,
  #                            maxIter=maxIter, verbose=verbose)
  #
  #   par.in <- c(as.vector(c(0,iter$alpha)),
  #               as.vector(c(0,beta.hat)),
  #               as.vector(iter$gamma),
  #               as.vector(iter$sdz),
  #               as.vector(sdy.hat))
  #
  #
  #
  #   # Number of covariates selected
  #   k <- length(which(iter$alpha!=0))
  #
  #   #HDIC <- 2*fn(par.in) + k*log(nX)
  #   #HDIC <- 2*fn(par.in) + k*pX^(1/3)
  #   HDIC <- 2*fn(par.in) + k*2*log(pX)
  #
  #
  #   return(HDIC)
  # }
  # parallel::stopCluster(cl)
  # options(warn = 0)
  #
  # HDIC <- cbind(lam0,HDIC.out)
  # bestlam <- HDIC[HDIC[,2]==min(HDIC[,2]),][1]

  #### Cross validation ####
  nfolds=10;
  foldid = sample(rep(seq(nfolds), length = nX))

  cl <- parallel::makeCluster(parallel::detectCores() - 1)
  doParallel::registerDoParallel(cl)
  cv.mse=foreach(i = 1:length(lam0), .combine=c, .export=c('stepsLassoSolver','stepsLeastR','stepsGLMNET','likelihoodSTEPS')) %dopar%{
    mse <- c()
    for(j in 1:nfolds){
      not.fold <- df3[foldid != j,]
      fold <- df3[foldid == j,]
      A.sub <- fold[,-c(1,2)]
      Y.sub <- fold[,1]
      Z.sub <- fold[,2]

      outlist = stepsLassoSolver(A=A.sub, Y2=Z.sub, Y1=Y.sub, X1=beta.hat, gamma=gam2, c1=c1, c2=c2,
                                 lambda = lam0[i], sigma2=sdz2, sigma1 = sdy.hat,
                                 method="BFGS", maxIter=maxIter, verbose=verbose)
      Ax=not.fold[,-c(1,2)] %*% outlist$alpha
      Axz=Ax-not.fold[,2];
      mse[j]=sum(Axz^2)
    }
    return(mean(mse))
  }
  parallel::stopCluster(cl)

  mse <- cbind(lam0,cv.mse)
  bestlam <- mse[mse[,2]==min(mse[,2]),][1]

  #### STEP 4 - Use stepsLassoSolver (stepsLeastR nested within) to iteratively estimate A ####
  options(warn = -1)
  estimates4 <- stepsLassoSolver(A=X, Y1=Y, X1=beta.hat, Y2=Z,
                                 c1=c1, c2=c2, lambda=bestlam, sigma1=sdy.hat, sigma2=sdz2, gamma=gam2,
                                 maxIter=maxIter, verbose=verbose)
  options(warn = 0)
  alpha4 <- estimates4$alpha[which(estimates4$alpha!=0),]

  if(length(alpha4)==0){
    return(list(beta.hat=beta.hat,
                sdy.hat=sdy.hat,
                initial.gamma=gam2,
                initial.sdz=sdz2,
                best.lambda=bestlam,
                alpha.lasso=NULL,
                X.sig.05=NULL,
                alpha.hat=NULL,
                sdz.hat=NULL,
                gamma.hat=NULL))
    invokeRestart("abort")
  }


  #### STEP 5 - Get D-Score p-value based on stepsLassoSolver estimates using stepsHDScoreTest() ####

  data.mat5 <- list(data=df3, c1=c1, c2=c2, sdy=sdy.hat, beta=beta.hat,
                    sdz=estimates4$sdz, gamma=estimates4$gamma, alpha=estimates4$alpha)

  alpha.HD.score <- stepsHDScoreTest(data.mat5)
  alpha5 <- alpha.HD.score[which(alpha.HD.score[,2]<0.05),2]
  S.hat <- noquote(names(alpha5))
  include5 <- which(alpha.HD.score[,2]<0.05)



  #### STEP 6 - Refit estimates with p<0.05 using stepsLD() ####

  data.mat6 <- list(c1=c1,c2=c2,data=cbind(Y,Z,X[,include5]))
  estimates6 <- stepsLD(data.mat6)
  rownames(estimates6$alpha.table) <- S.hat
  colnames(estimates6$alpha.table)[1] <- "alpha.hat"

  #### STEP 7 - Return estimates ####


  list(beta.hat=beta.hat,
       sdy.hat=sdy.hat,
       initial.gamma=gam2,
       initial.sdz=sdz2,
       best.lambda=bestlam,
       alpha.lasso=alpha.HD.score,
       X.sig.05=S.hat,
       alpha.hat.table=estimates6$alpha.table,
       alpha.hat=estimates6$alpha.hat,
       sdz.hat=estimates6$sdz.hat,
       gamma.hat=estimates6$gamma.hat)



}






