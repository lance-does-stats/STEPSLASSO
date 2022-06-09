#' Perform Decorrelated Score test based on Ning et al. for ST-EPS data
#'
#' @param data.mat List containing Y, Z, X, c1, c2. Required.
#' @param position indicates position in coefficient matrix of variable to test. Required.
#' @param null_value the value to test the null hypothesis against (default=0)
#' @import magrittr glmnet
#' @export
#' @example


stepsHDScoreTest <- function(list, null_value=0, TestParallel=F, reserveNcores=1){

  c1 <- list$c1
  c2 <- list$c2
  df <- list$data

  Y <- df[,1]
  X <- df[,-c(1,2)]
  Z <- df[,2]

  nX <- dim(X)[1]
  pX <- dim(X)[2]

  c1=list$c1
  c2=list$c2
  theta <- sqrt(list$gamma^2 * list$sdz^2 + list$sdy^2) %>% as.numeric()
  #A0 <- estimates$alpha0
  A <- list$alpha
  #B0 <- 0
  B <- list$beta
  gam <- list$gamma
  sdz <- list$sdz
  sdy <- list$sdy

  pval_score=c()
  #alpha.correct=c()
  if (TestParallel==T){
    cl <- parallel::makeCluster(parallel::detectCores() - reserveNcores)
    doParallel::registerDoParallel(cl)
    pval_score = foreach(i=1:pX, .combine=c, .export=c('cv.glmnet')) %dopar%{

      As <- A; As[i] <- null_value
      if(length(A)==1){
        n <-  length(X)
        K <-  1
        XA <-  X*As
        XB <-  X*B
      }else{
        n <-  nrow(X)
        K <-  ncol(X)
        XA <- X%*%As
        XB <- X%*%B
      }

      t1 <- c1 - (XB + gam*XA)# lower
      t2 <- c2 - (XB + gam*XA)# upper

      f1 <- dnorm(t1, sd=theta)
      f2 <- dnorm(t2, sd=theta)
      F1 <- pnorm(t1, sd=theta)
      F2 <- pnorm(t2, sd=theta)
      P.S1 <- F1 + 1 - F2

      mv=(f1-f2)/P.S1
      vv=(t1*f1-t2*f2)/P.S1 + mv^2
      # Score function
      SM <- (t(X)%*%(((Z-XA)/sdz^2) + (f1 - f2)/P.S1*gam/theta))/nX

      # Fisher's information matrix
      J = matrix(0,nrow=K,ncol=K)
      J[1:K, 1:K]<-t(X) %*%diag(as.numeric(-1/sdz^2 + (gam^2 / theta^2)*(vv)))%*%X/nX
      #Jinv = solve(J)

      HM <- -J

      #### Fit w with Lasso ####

      W <- matrix(0, pX, pX-1);
      alpha_i <- As


      la  = rep(0,nX)                           # Gradient w.r.t parameter of interest
      lb  = matrix(0,nrow = nX, ncol=pX-1)
      for(k in 1:nX){
        la[k]=(Z[k]-X[k,]%*%alpha_i)*X[k,i]/sdz^2 - as.numeric(mv[k]*X[k,i]*gam/theta)
        lb[k,]=as.vector(Z[k]-X[k,]%*%alpha_i)*X[k,-i]/sdz^2 - as.numeric(mv[k]*X[k,-i]*gam/theta)

      }

      cv.fit = cv.glmnet(lb,la,standardize = FALSE,intercept = FALSE)
      fit     = cv.fit$glmnet.fit
      tmp     = which(fit$lambda == cv.fit$lambda.min)

      if (sum(fit$beta[,tmp]!=0)==0){
        W[i,] = rep(0,pX-1)
      } else {
        W[i,]    = fit$beta[,tmp]
      }

      if (i == 1){
        var   = max(HM[i,i] - W[i,]%*%HM[c((i+1):pX),i])
      } else if (i == pX){
        var   = max(HM[i,i] - W[i,]%*%HM[c(1:(i-1)),i])
      } else {
        var   = max(HM[i,i] - W[i,]%*%HM[c(1:(i-1),(i+1):pX),i])
      }

      ## Decorrelated Score
      S = SM[i,] - W[i,] %*% SM[-i,]
      #pval_score[i]=2*pnorm(-abs(sqrt(nX)*S/sqrt(max(var,0.1))))

      pval_s=1 - pchisq(n*S^2/var, df=1)
      #alpha.correct[i] <- A[i] - S/var
      return(pval_s)
    } #dopar end

    parallel::stopCluster(cl)

    } else{

      for(i in 1:pX){
        As <- A; As[i] <- null_value

        if(length(A)==1){
          n <-  length(X)
          K <-  1
          XA <-  X*As
          XB <-  X*B
        }else{
          n <-  nrow(X)
          K <-  ncol(X)
          XA <- X%*%As
          XB <- X%*%B
        }

        t1 <- c1 - (XB + gam*XA)# lower
        t2 <- c2 - (XB + gam*XA)# upper

        f1 <- dnorm(t1, sd=theta)
        f2 <- dnorm(t2, sd=theta)
        F1 <- pnorm(t1, sd=theta)
        F2 <- pnorm(t2, sd=theta)
        P.S1 <- F1 + 1 - F2

        mv=(f1-f2)/P.S1
        vv=(t1*f1-t2*f2)/P.S1 + mv^2
        # Score function
        SM <- (t(X)%*%(((Z-XA)/sdz^2) + (f1 - f2)/P.S1*gam/theta))/nX

        # Fisher's information matrix
        J = matrix(0,nrow=K,ncol=K)
        J[1:K, 1:K]<-t(X) %*%diag(as.numeric(-1/sdz^2 + (gam^2 / theta^2)*(vv)))%*%X/nX
        #Jinv = solve(J)

        HM <- -J

        #### Fit w with Lasso ####

        W <- matrix(0, pX, pX-1);
        alpha_i <- As


        la  = rep(0,nX)                           # Gradient w.r.t parameter of interest
        lb  = matrix(0,nrow = nX, ncol=pX-1)
        for(k in 1:nX){
          la[k]=(Z[k]-X[k,]%*%alpha_i)*X[k,i]/sdz^2 - as.numeric(mv[k]*X[k,i]*gam/theta)
          lb[k,]=as.vector(Z[k]-X[k,]%*%alpha_i)*X[k,-i]/sdz^2 - as.numeric(mv[k]*X[k,-i]*gam/theta)

        }

        cv.fit = cv.glmnet(lb,la,standardize = FALSE,intercept = FALSE)
        fit     = cv.fit$glmnet.fit
        tmp     = which(fit$lambda == cv.fit$lambda.min)

        if (sum(fit$beta[,tmp]!=0)==0){
          W[i,] = rep(0,pX-1)
        } else {
          W[i,]    = fit$beta[,tmp]
        }

        if (i == 1){
          var   = max(HM[i,i] - W[i,]%*%HM[c((i+1):pX),i])
        } else if (i == pX){
          var   = max(HM[i,i] - W[i,]%*%HM[c(1:(i-1)),i])
        } else {
          var   = max(HM[i,i] - W[i,]%*%HM[c(1:(i-1),(i+1):pX),i])
        }

        ## Decorrelated Score
        S = SM[i,] - W[i,] %*% SM[-i,]
        #pval_score[i]=2*pnorm(-abs(sqrt(nX)*S/sqrt(max(var,0.1))))

        pval_score[i]=1 - pchisq(n*S^2/var, df=1)
      } #for end
  } #else end

  return(cbind(A, pval_score))
} #function end



