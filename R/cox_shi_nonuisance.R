#' @description Run the Cox and Shi test for the null E[Y] <= 0 for Y~N(mu,Sigma)
#'
#'
cox_shi_nonuisance <- function(Y, sigma, alpha = 0.05){
  #We want to minimize (Y-mu)' sigma^{-1} (Y-mu) s.t mu <=0
  # This is the same as minimizing mu' Sigma^-1 mu - 2 (Sigma^-1 Y)' mu + Y'Sigma^{-1} Y
  sigmaInv <- solve(sigma)

  #Convert Y to a column vector if not already
  Y <- as.matrix(Y,ncol = length(Y))

  #Quadprog solves min 1/2 x'Dx - d'x s.t. A'x >= b
  # Thus, we set D = 2 * sigmaInv,
  # d= 2*Sigma^{-1} Y + Sigma^{-1} Y', A=-I, b=0
  # Then we add Y'Sigma^{-1} Y to the minimum

  qp <-
  quadprog::solve.QP(Dmat = 2*sigmaInv,
                     dvec = 2 * sigmaInv %*% Y,
                     Amat = -diag(NROW(sigmaInv)),
                     bvec = matrix(0,nrow = NROW(sigmaInv), ncol =1))

  test_stat <- qp$value + t(Y) %*% sigmaInv %*% Y

  #Find which constraints are binding, up to a tolerance
  binding_index <- which(abs(qp$solution)<10^-5)

  chisquared_df <- length(binding_index)

  if(chisquared_df == 0){
    return(list(reject = 0,
                test_stat = test_stat,
                cv = 0))
  }else{
    cv <- qchisq(p = 1-alpha,
                 df = length(binding_index))

    return(list(reject = test_stat > cv,
                test_stat = test_stat,
                cv = cv))

  }



}
