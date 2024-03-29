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


stepsLasso <- function(Y, c1, c2, Z, X, beta.hat, sdy.hat, method=NULL,lam0=NULL, maxIter=1000, verbose=FALSE, LassoParallel=T, TestParallel=F, reserveNcores=1){
  #### STEP 1 - Estimate Beta.hat and sdy.hat from EPS-LASSO ####

  beta.hat <- beta.hat
  sdy.hat <- sdy.hat



  #### STEP 2 - Get initial estimate for sdz and gamma from GLMNET ####

  data.mat2 <- list(c1=c1, c2=c2, data=cbind(Y,Z,X), beta=beta.hat, sd.y=sdy.hat)
  estimates2 <- stepsGLMNET(data.mat2)


  # initial estimate for sdz
  sdz2 <- estimates2$int.sdz
  gam2 <- estimates2$int.gamma
  optim2Worked <- estimates2$optimized


  #### STEP 3 - Choose best lambda to pass to stepsLassoSolver ####

  nX <-  nrow(X)
  pX <-  ncol(X)


  # Create range of lambda for optimization
  lam0=lam0
  if(is.null(lam0)){
    lamMax=max(abs(t(X)%*%Z))
    if(nX<=pX){
      lamMin=0.01*lamMax
    }else{
      lamMin=0.0001*lamMax
    }
    lamR=exp((log(lamMax)-log(lamMin))/99)
    lam0=c(lamMax,lamMax*lamR^(-seq(1:99)))
    lam0=lam0[order(lam0)];
  }

  df3 <- cbind(Y,Z,X) %>% as.matrix()
  #### HDIC or CV for bestlam ####
  if(method=="HDIC"){
    ### HDIC ####
    K=pX
    fn = function(par){
      para.ll = list(A0= par[1], A = par[1+1:K],
                     B0=par[2+K], B = par[2+K+1:K],
                     gam = par[3+2*K],
                     sd.z = par[4+2*K],
                     sd.y = par[5+2*K])
      res.ll = likelihoodSTEPS(df3, c1=c1, c2=c2, para.ll)
      return(res.ll$ll)
    }

    # Parallel
    options(warn = -1)
    cl <- parallel::makeCluster(parallel::detectCores() - reserveNcores)
    doParallel::registerDoParallel(cl)
    HDIC.out=foreach(i = 1:length(lam0), .combine=c, .export=c('stepsLassoSolver','stepsLeastR','stepsGLMNET','likelihoodSTEPS'), .packages=c('MASS','glmnet')) %dopar%{
      HDIC <- c()


      iter <- stepsLassoSolver(A=X, Y2=Z, Y1=Y, X1=beta.hat, gamma=gam2, c1=c1, c2=c2,
                               lambda = lam0[i], sigma2=sdz2, sigma1 = sdy.hat,
                               maxIter=maxIter, verbose=verbose)

      par.in <- c(as.vector(c(0,iter$alpha)),
                  as.vector(c(0,beta.hat)),
                  as.vector(iter$gamma),
                  as.vector(iter$sdz),
                  as.vector(sdy.hat))



      # Number of covariates selected
      k <- length(which(iter$alpha!=0))

      #HDIC <- 2*fn(par.in) + k*log(nX)
      HDIC <- 2*fn(par.in) + k*pX^(1/3)
      #HDIC <- 2*fn(par.in) + k*2*log(pX)


      return(HDIC)
    }
    parallel::stopCluster(cl)
    options(warn = 0)

    HDIC <- cbind(lam0,HDIC.out)
    bestlam <- HDIC[HDIC[,2]==min(HDIC[,2]),][1]
    method.out="HDIC"
  } else{

    ### Cross validation ####
    if(nX<100){
      nfolds=3
    } else{
      nfolds=10;
    }

    foldid = sample(rep(seq(nfolds), length = nX))
    if(LassoParallel==T){
      cl <- parallel::makeCluster(parallel::detectCores() - reserveNcores)
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
                                     lambda = lam0[i], sigma2=sdz2, sigma1 = sdy.hat, maxIter=maxIter, verbose=verbose)
          Ax=not.fold[,-c(1,2)] %*% outlist$alpha
          Axz=Ax-not.fold[,2];
          mse[j]=sum(Axz^2)
        }
        return(mean(mse))
      }
      parallel::stopCluster(cl)

    } else{

      for(i in 1:length(lam0)){
        mse=c()
        for(j in 1:nfolds){
          not.fold <- df3[foldid != j,]
          fold <- df3[foldid == j,]
          A.sub <- fold[,-c(1,2)]
          Y.sub <- fold[,1]
          Z.sub <- fold[,2]

          outlist = stepsLassoSolver(A=A.sub, Y2=Z.sub, Y1=Y.sub, X1=beta.hat, gamma=gam2, c1=c1, c2=c2,
                                     lambda = lam0[i], sigma2=sdz2, sigma1 = sdy.hat, maxIter=maxIter, verbose=verbose)
          Ax=not.fold[,-c(1,2)] %*% outlist$alpha
          Axz=Ax-not.fold[,2];
          mse[j]=sum(Axz^2)
        }
      }
    }


    mse <- cbind(lam0,cv.mse)
    bestlam <- mse[mse[,2]==min(mse[,2]),][1]
    method.out="CV"
  }




  #### STEP 4 - Use stepsLassoSolver (stepsLeastR nested within) to iteratively estimate A ####
  options(warn = -1)
  estimates4 <- stepsLassoSolver(A=X, Y1=Y, X1=beta.hat, Y2=Z,
                                 c1=c1, c2=c2, lambda=bestlam, sigma1=sdy.hat, sigma2=sdz2, gamma=gam2,
                                 maxIter=maxIter, verbose=verbose)
  options(warn = 0)
  alpha4 <- estimates4$alpha[which(estimates4$alpha!=0),]

  if(length(alpha4)==0){
      list(beta.hat=beta.hat,
                  sdy.hat=sdy.hat,
                  initial.gamma=gam2,
                  initial.sdz=sdz2,
                  best.lambda=bestlam,
                  alpha.refit=NULL,
                  X.selected=NULL,
                  sdz.hat=NULL,
                  gamma.hat=NULL,
                  optim2Worked=optim2Worked,
                  optim6Worked=1,
                  method.out=method.out)
  } else{
    S.hat <- noquote(names(alpha4))

      #### STEP 7 - Return estimates ####
      list(beta.hat=beta.hat,
           sdy.hat=sdy.hat,
           initial.gamma=gam2,
           initial.sdz=sdz2,
           best.lambda=bestlam,
           alpha.refit=NULL,
           X.selected=S.hat,
           sdz.hat=sdz2,
           gamma.hat=gam2,
           optim2Worked=optim2Worked,
           optim6Worked=1,
           method.out=method.out)
  }
  # if(length(alpha4)==0){
  #   #### STEP 5 - Return null values since no X's selected ####
  #   A <- rep(0,pX)
  #   pval_score <- rep(1,pX)
  #   list(beta.hat=beta.hat,
  #               sdy.hat=sdy.hat,
  #               initial.gamma=gam2,
  #               initial.sdz=sdz2,
  #               best.lambda=bestlam,
  #               alpha.refit=cbind(A,pval_score),
  #               X.selected=NULL,
  #               sdz.hat=NULL,
  #               gamma.hat=NULL,
  #               optim2Worked=optim2Worked,
  #               optim6Worked=1)
  #
  # } else if (length(alpha4)==1 || length(alpha4)==2){
  #   S.hat <- noquote(names(alpha4))
  #   include5 <- which(estimates4$alpha!=0) + 2
  #
  #
  #   #### STEP 5 - Refit estimates with p<0.05 using stepsLD() ####
  #   data.mat5 <- list(data=df3[,c(1,2,include5)], c1=c1, c2=c2, sd.y=sdy.hat, B=beta.hat[include5-2])
  #   estimates5 <- optimSTEPS(data.mat5)
  #
  #   #### STEP 6 - Get LD-Score p-value based on scoreTestSTEPS ####
  #   alpha5 <- estimates5$alpha
  #   names(alpha5) <- S.hat
  #
  #   if(length(alpha5)==1){
  #     pval_score <- scoreTestSTEPS(data.mat5)$`Score P-value`
  #   } else{
  #     pval_score <- c(scoreTestSTEPS(data.mat5, position = 1)$`Score P-value`,
  #                     scoreTestSTEPS(data.mat5, position = 2)$`Score P-value`)
  #   }
  #   alpha.LD.score=cbind(A=alpha5,pval_score)
  #
  #
  #   #### STEP 7 - Return estimates ####
  #   list(beta.hat=beta.hat,
  #        sdy.hat=sdy.hat,
  #        initial.gamma=gam2,
  #        initial.sdz=sdz2,
  #        best.lambda=bestlam,
  #        alpha.refit=alpha.LD.score,
  #        X.selected=S.hat,
  #        sdz.hat=estimates5$sdz,
  #        gamma.hat=estimates5$gamma,
  #        optim2Worked=optim2Worked,
  #        optim6Worked=1)
  #
  # } else if(length(alpha4)>(pX-3)){
  #
  #   #### STEP 5 -  Use alpha estimates from stepsLassoSolver and bestlam####
  #   S.hat <- noquote(names(alpha4))
  #   include5 <- which(estimates4$alpha!=0) + 2
  #   alpha5 <- alpha4
  #   df5 <- df3[,c(1,2,include5)]
  #   beta.selected <- beta.hat[include5-2]
  #   data.mat5 <- list(data=df5, c1=c1, c2=c2, sd.y=sdy.hat, B=beta.selected)
  #
  #   #### STEP 6 -  Get LD-Score p-value based on stepsLD() ####
  #   uni.pval=matrix(nrow= length(alpha5),ncol = 2)
  #   for(j in 1:length(alpha5)){
  #     data.mat.uni <- data.mat5
  #     data.mat.uni$data <- data.mat.uni$data[,c(1,2,j+2)]
  #     uni.pval[j,] <- stepsLD(data.mat.uni)$alpha.table[,1:2]
  #   }
  #
  #   rownames(uni.pval) <- paste0('X',1:length(alpha5))
  #   colnames(uni.pval) <- c('A','pval_score')
  #
  #
  #   #### STEP 7 - Return estimates ####
  #   list(beta.hat=beta.hat,
  #        sdy.hat=sdy.hat,
  #        initial.gamma=gam2,
  #        initial.sdz=sdz2,
  #        best.lambda=bestlam,
  #        alpha.refit=uni.pval,
  #        X.selected=S.hat,
  #        sdz.hat=estimates4$sdz,
  #        gamma.hat=estimates4$gamma,
  #        optim2Worked=optim2Worked,
  #        optim6Worked=0)
  #
  #
  #
  # } else{
  #   S.hat <- noquote(names(alpha4))
  #   include5 <- which(estimates4$alpha!=0) + 2
  #
  #
  #   #### STEP 5 -  Refit estimates with p<0.05 using stepsLD() ####
  #   df5 <- df3[,c(1,2,include5)]
  #   beta.selected <- beta.hat[include5-2]
  #   data.mat5 <- list(data=df5, c1=c1, c2=c2, sd.y=sdy.hat, B=beta.selected)
  #
  #   estimates5 <- optimSTEPS(data.mat5)
  #   alpha5 <- estimates5$alpha
  #   names(alpha5) <- S.hat
  #
  #
  #   #### STEP 6 -  Get HD-Score p-value based on stepsLassoSolver estimates using stepsHDScoreTest()####
  #   data.mat6 <- list(data=df5, c1=c1, c2=c2, sdy=sdy.hat, beta=beta.selected,
  #                     sdz=estimates5$sdz, gamma=estimates5$gamma, alpha=alpha5)
  #
  #   alpha.HD.score <- stepsHDScoreTest(data.mat6, TestParallel=TestParallel, reserveNcores=reserveNcores)
  #
  #
  #   #### STEP 7 - Return estimates ####
  #   list(beta.hat=beta.hat,
  #        sdy.hat=sdy.hat,
  #        initial.gamma=gam2,
  #        initial.sdz=sdz2,
  #        best.lambda=bestlam,
  #        alpha.refit=alpha.HD.score,
  #        X.selected=S.hat,
  #        sdz.hat=estimates5$sdz,
  #        gamma.hat=estimates5$gamma,
  #        optim2Worked=optim2Worked,
  #        optim6Worked=1)
  # }


}






