#' This function solves EPSLASSO with an initial guess of beta and sigma.
#'
#' @param A Matrix of predictors. Required.
#' @param y1 Response vector of primary trait. Required.
#' @param x1 Vector of coefficients for primary trait. Required.
#' @param y2 Response vector of secondary trait. Required.
#' @param z L_1 norm regularization parameter (z >=0). Required.
#' @param c1 Left censored point. Required.
#' @param c2 Right censored point. Required.
#' @param sigma1 Estimate of sdy (PT) Required.
#' @param sigma2 Estimate of sdz (ST) Required.
#' @param gamma Estimate of gamma eps(Y)~Z. Required.
#' @param maxIter Maximum iteration number for approximation procedure. Default is 1000.
#' @param verbose Print debugging info or not.
#' @import
#' @export
#' @example
#' # data.mat <- generateSTEPS(n=1000, p=100, q.a=3, q.b=3, A_coef=0.7, B_coef=0.7, overlap=0,
#'                          gam=0.8, sd.y=1, sd.z=1, eps=0.2)
#'
#' int.estimates <- stepsGLMNET(data.mat)
#' int.sigma2 <- as.numeric(int.estimates$int.sdz)
#' int.gamma <- as.numeric(int.estimates$int.gamma)
#' B.hat <- data.mat$B
#' sdy.hat <- data.mat$sd.y
#
#' temp2 <- stepsLassoSolver(A=data.mat$data[,-c(1,2)],
#'                      Y1=data.mat$data[,1], Y2=data.mat$data[,2], X1=B.hat,
#'                      lambda=6,
#'                      c1=data.mat$c1, c2=data.mat$c2,
#'                      sigma1=sdy.hat, sigma2=int.sigma2, gamma=int.gamma)


stepsLassoSolver=function(A, Y1, X1, Y2, c1, c2, lambda, sigma1, sigma2, gamma, method="L-BFGS-B", maxIter=1000, verbose=FALSE){

  pX=ncol(A)
  nX=nrow(A)
  x=matrix(0,pX,1);
  sigma1=as.numeric(sigma1)
  sigma2=as.numeric(sigma2)
  gamma=as.numeric(gamma)

  esp_beta=1
  esp_var=1
  esp_gam=1

  fit=stepsLeastR(A,y1=Y1, y2=Y2, x1=X1, z=lambda,
                  c1=c1, c2=c2,
                  sigma1=sigma1, sigma2=sigma2, gamma=gamma)
  x=fit$x
  step=0

  while ( (esp_beta>0.001 || esp_var>0.001 || esp_gam>0.001) & step<=maxIter ){
    step=step+1
    sigma_old=sigma2
    x_old=x
    gamma_old=gamma

    ##### ESTIMATE SIGMA2 using optim() function ####
    par.in <- c(as.vector(c(0,x)),
                as.vector(c(0,X1)),
                as.vector(gamma),
                as.vector(sigma2),
                as.vector(sigma1))

    K <- pX
    df <- cbind(Y1,Y2,A)
    fn = function(par){
      para.ll = list(A0= par[1], A = par[1+1:pX],
                     B0=par[2+pX], B = par[2+pX+1:pX],
                     gam = par[3+2*pX],
                     sd.z = par[4+2*pX],
                     sd.y = par[5+2*pX])
      res.ll = likelihoodSTEPS(df, c1=c1, c2=c2, para.ll)
      return(res.ll$ll)
    }

    gr = function(par){
      para.ll = list(A0= par[1], A = par[1+1:pX],
                     B0=par[2+pX], B = par[2+pX+1:pX],
                     gam = par[3+2*pX],
                     sd.z = par[4+2*pX],
                     sd.y = par[5+2*pX])
      res.ll = likelihoodSTEPS(df, c1=c1, c2=c2, para.ll)
      return(res.ll$d.ll)
    }

    if(method=="CG" || method == "BFGS"){
      res.opt=tryCatch(optim(par.in,fn,gr,method = method), error=function(e) optim(par.in,fn,gr,method = "L-BFGS-B"))
    } else{
      res.opt=tryCatch(optim(par.in, fn, gr, method = method,
                    lower=c(0, rep(0,pX), 0, rep(0,pX), 0, 0.0001, 0.0001),
                    upper=c(10000, rep(10000,pX), 10000, rep(10000,pX), 1000, 1000, 1000)), optim(par.in,fn,gr,method = "CG"))
    }

    sigma2=res.opt$par[4+2*pX]
    gamma=res.opt$par[3+2*pX]

    #### If no change, then stop
    if(sigma2==sigma_old || gamma==gamma_old){
      break;
    } ## otherwise, refit alphas based on new sigma2
    fit=stepsLeastR(A,y1=Y1, y2=Y2, x1=X1,z=lambda,
                    c1=c1, c2=c2,
                    sigma1=sigma1, sigma2=sigma2, gamma=gamma)
    if(fit$ill_tag==1){
      sigma2=sigma_old
      x=x_old
      gamma=gamma_old
      break
    }else{
      x=fit$x
    }

    esp_beta=norm(x-x_old)
    esp_var=abs(sigma2-sigma_old)
    esp_gam=abs(gamma-gamma_old)
    if(step>1 && esp_beta>pre_esp_beta && esp_var>pre_esp_var && esp_gam>pre_esp_gam){
      if(verbose){
        print("Diverging.")
      }
      sigma2=sigma_old
      x=x_old
      gamma=gamma_old
      break;
    }
    pre_esp_beta = esp_beta
    pre_esp_var = esp_var
    pre_esp_gam = esp_gam
  }
  out = list(lambda = lambda, sdz=sigma2, gamma=gamma, alpha=x)
}

