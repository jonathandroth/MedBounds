#' @description Run the Cox and Shi test for the null E[Y] <= 0 for Y~N(mu,Sigma)
#'
#'
cox_shi_nonuisance <- function(Y, sigma, alpha = 0.05){

  if(min(base::eigen(sigma, only.values = T)$values) < 10^-6){
    #If sigma is not full-rank we extract the full rank component

    eigendecomp <- eigen(sigma)
    eigenvals <- eigendecomp$values
    eigenvecs <- eigendecomp$vectors

    tol <- 10^-6
    positive_indices <- which(eigenvals > tol)

    #Create a matrix that selects rows with positive indices
    M <- diag(length(eigenvals))
    M <- as.matrix(M[positive_indices,])

    #Deal with annoying fact that R converts to column vector if lenght 1
    if(length(positive_indices) == 1){M <- t(M)}

    Xstar <- M %*% t(eigenvecs) %*% Y
    V_Xstar <- diag( x = eigenvals[positive_indices],
                     nrow = length(eigenvals[positive_indices]))


    #Matrix for converting between Xstar and Y
    A <- eigenvecs %*% t(M)

    #Find constant difference between Y and A Xstar
    #This becomes the constant in your constraints, i.e. A mu <= b
    b <- A %*% Xstar - Y
    Y <- Xstar
    sigma <- V_Xstar
  }else{
    A <- diag(NROW(sigma))
    b <- matrix(0,nrow = NROW(A), ncol =1)
  }
  #We want to minimize (Y-mu)' sigma^{-1} (Y-mu) s.t mu <=0
  # This is the same as minimizing mu' Sigma^-1 mu - 2 (Sigma^-1 Y)' mu + Y'Sigma^{-1} Y
  sigmaInv <- solve(sigma)

  #Convert Y to a column vector if not already
  Y <- as.matrix(Y,ncol = length(Y))

  #Quadprog solves min 1/2 x'Dx - d'x s.t. A'x >= b
  # Thus, we set D = 2 * sigmaInv,
  # d= 2*Sigma^{-1} Y + Sigma^{-1} Y', b=0
  # Then we add Y'Sigma^{-1} Y to the minimum

  qp <-
  quadprog::solve.QP(Dmat = 2*sigmaInv,
                     dvec = 2 * sigmaInv %*% Y,
                     Amat = -t(A),
                     bvec = -b)

  test_stat <- qp$value + t(Y) %*% sigmaInv %*% Y

  #Find which constraints are binding, up to a tolerance
  binding_index <- which(abs(A %*% qp$solution - b)<10^-5)

  chisquared_df <- length(binding_index)

  if(chisquared_df == 0){
    return(list(reject = 0,
                test_stat = test_stat,
                cv = 0,
                pval = 1))
  }else{
    cv <- qchisq(p = 1-alpha,
                 df = length(binding_index))

    pval <- pchisq(q = test_stat,
                   lower.tail = FALSE,
                   df = length(binding_index))

    return(list(reject = test_stat > cv,
                test_stat = test_stat,
                cv = cv,
                pval = pval))

  }



}
