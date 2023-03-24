#' Test for high-dimensional sparse regression for secondary trait analysis under extreme phenotype sampling data.
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



stepsLasso <- function(Y, c1, c2, Z, X, B0.hat, beta.hat, sdy.hat, selection="CV10", method="BFGS", maxIter=1000, verbose=FALSE){

  #### STEP 1 - Estimate Beta.hat and sdy.hat from EPS-LASSO ####
  B0 <- B0.hat
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
    lamMax=max(abs(t(X)%*%Y))
    if(nX<pX){
      lamMin=0.01*lamMax
    }else{
      lamMin=0.0001*lamMax
    }
    lamR=exp((log(lamMax)-log(lamMin))/99)
    lam0=c(lamMax,lamMax*lamR^(-seq(1:99)))
    lam0=lam0[order(lam0)];
  }

  # Function to calculate -logLikelihood value
  df <- cbind(Y,Z,X)

  if(selection %in% c("AIC", "BIC", "HDIC1","HDIC2", "EBIC","HQC")){
    #  #### HDICs ####
    fn = function(par){
      para.ll = list(A0= par[1], A = par[1+1:pX],
                     B0=par[2+pX], B = par[2+pX+1:pX],
                     gam = par[3+2*pX],
                     sd.z = par[4+2*pX],
                     sd.y = par[5+2*pX])
      res.ll = likelihoodSTEPS(df, c1=c1, c2=c2, para.ll)
      return(res.ll$ll.sum)
    }

    # Parallel
    options(warn = -1)
    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    HDIC.out=foreach(i = 1:length(lam0), .combine=c, .export=c('stepsLassoSolver','stepsLeastR','stepsGLMNET','likelihoodSTEPS'), .packages=c('MASS','glmnet')) %dopar%{
      HDIC <- c()

      intercept <- rep(1,nX)
      X.int <- cbind(intercept,X)
      beta.hat.int <- c(B0, beta.hat)

      iter <- stepsLassoSolver(A=X.int, Y2=Z, Y1=Y, X1=beta.hat.int, gamma=gam2, c1=c1, c2=c2,
                               lambda = lam0[i], sigma2=sdz2, sigma1 = sdy.hat, method=method,
                               maxIter=maxIter, verbose=verbose)

      par.in <- c(as.vector(c(iter$alpha)),
                  as.vector(c(beta.hat.int)),
                  as.vector(iter$gamma),
                  as.vector(iter$sdz),
                  as.vector(sdy.hat))

      # Number of covariates selected
      k <- length(which(iter$alpha!=0))

      if(selection == "AIC"){
        HDIC <- 2*fn(par.in) + 2*k

      } else if(selection == "BIC") {
        HDIC <- 2*fn(par.in) + k*log(nX)

      } else if(selection == "HDIC1") {
        HDIC <- 2*fn(par.in) + k*pX^(1/3)

      } else if(selection == "HDIC2") {
        HDIC <- 2*fn(par.in) + k*2*log(pX)

      } else if(selection == "EBIC") {
        HDIC <- 2*fn(par.in) + k*log(nX) + 0.1*k*log(pX)

      } else if(selection == "HQC") {
        HDIC <- 2*fn(par.in) + log(log(nX))
      }


      return(HDIC)
    }
    parallel::stopCluster(cl)
    options(warn = 0)

    HDIC <- cbind(lam0,HDIC.out)
    #plot(HDIC[,1], HDIC[,2])
    bestlam <- HDIC[HDIC[,2]==min(HDIC[,2]),][1]

  } else{

    #### Cross validation ####
    if(selection == "CV10"){
      nfolds=10;
    } else if (selection == "CV5"){
      nfolds=5;
    }

    foldid = sample(rep(seq(nfolds), length = nX))

    cl <- parallel::makeCluster(parallel::detectCores() - 1)
    doParallel::registerDoParallel(cl)
    cv.mse=foreach(i = 1:length(lam0), .combine=c, .export=c('stepsLassoSolver','stepsLeastR','stepsGLMNET','likelihoodSTEPS')) %dopar%{
      mse <- c()
      for(j in 1:nfolds){
        not.fold <- df[foldid != j,]
        fold <- df[foldid == j,]
        A.sub <- fold[,-c(1,2)]
        Y.sub <- fold[,1]
        Z.sub <- fold[,2]

        outlist = stepsLassoSolver(A=A.sub, Y2=Z.sub, Y1=Y.sub, X1=beta.hat, gamma=gam2, c1=c1, c2=c2,
                                   lambda = lam0[i], sigma2=sdz2, sigma1 = sdy.hat,
                                   maxIter=maxIter, verbose=verbose)
        Ax=not.fold[,-c(1,2)] %*% outlist$alpha
        Axz=Ax-not.fold[,2];
        mse[j]=sum(Axz^2)
      }
      return(mean(mse))
    }
    parallel::stopCluster(cl)

    mse <- cbind(lam0,cv.mse)
    #plot(mse[,1], mse[,2])
    bestlam <- mse[mse[,2]==min(mse[,2]),][1]
  }




  #### STEP 4 - Use stepsLassoSolver (stepsLeastR nested within) to iteratively estimate A ####
  #method=method.in <- "CG"
  estimates4 <- stepsLassoSolver(A=X, Y1=Y, X1=beta.hat, Y2=Z,
                                 c1=c1, c2=c2, lambda=bestlam, sigma1=sdy.hat, sigma2=sdz2, gamma=gam2,
                                 method=method, maxIter=maxIter, verbose=verbose)

  alpha4 <- estimates4$alpha[which(estimates4$alpha!=0),]
  alpha4

  if(length(alpha4)==0){
    return(  list(included.X=0,
                  included.X.index= NULL,
                  sdz.hat=estimates4$sdz,
                  gamma.hat=estimates4$gamma,
                  sdy.hat=sdy.hat,
                  alpha.hat=alpha4,
                  beta.hat=beta.hat,
                  best.lambda=bestlam,
                  num.p.selected= 0,
                  num.p.entered= pX,
                  prop.p.reduced= 0,
                  selection=selection)
    )
    invokeRestart("abort")
  }
  S.hat <- noquote(names(alpha4))
  include4 <- which(estimates4$alpha!=0)


  #### STEP 5 - Return estimates ####

  list(included.X=S.hat,
       included.X.index= include4,
       sdz.hat=estimates4$sdz,
       gamma.hat=estimates4$gamma,
       sdy.hat=sdy.hat,
       alpha.hat=alpha4,
       beta.hat=beta.hat,
       best.lambda=bestlam,
       num.p.selected= length(include4),
       num.p.entered= pX,
       prop.p.reduced=1 - length(include4)/pX,
       selection=selection)



}









