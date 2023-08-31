#' @title Hypothesis test for the sharp null of full mediation
#' @description This function tests the sharp null that Y(1,m) = Y(0,m). The
#'   outcome and mediator are both assumed to take finitely many different
#'   values. Several options are available for infernece.
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable, which is assumed to take a discrete
#'   support
#'   @param method Method to use. One of "ARP, "CS", "FSST", "CR"
#' @param ordering A list with length equal to the cardinality of the support of
#'   the mediator variable. The name of each element corresponds to a point in
#'   the support, and each element is a vector that collects all m values that
#'   are less than or equal to this point. If ordering = NULL, the standard
#'   ordering is used.
#' @param B Bootstrap size, default is 500
#' @param cluster Cluster for bootstrap
#' @param weight.matrix Weight matrix used to implement FSST. Possible options
#'   are "diag", "avar", "identity." Defaults is "diag" as in FSST.
#' @param hybrid_kappa The first-stage size value of the ARP hybrid test. If NULL, the ARP conditional (non-hybrid) test is used. Default is alpha/10
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
#' @param alpha Significance level. Default is 0.05
#' @export
#'

test_sharp_null <- function(df,
                            d,
                            m,
                            y,
                            method,
                            ordering = NULL,
                            B = 500,
                            cluster = NULL,
                            weight.matrix = "diag",
                            hybrid_kappa = alpha/10,
                            num_Ybins = NULL,
                            alpha = 0.05){

  ## Process the inputted df ----

  #Remove missing values

  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)

  #Discretize y if needed
  if(!is.null(num_Ybins)){
    df[[y]] <- discretize_y(yvec = df[[y]], numBins = num_Ybins)
  }

  yvec <- df[[y]]
  dvec <- df[[d]]
  mvec <- df[[m]]

  ## Construct the A matrices and beta.shp ----

  #Specify whether method requires us to input only inequalities
  inequalities_only <- ifelse(method %in% c("ARP","CS"),
                              TRUE, FALSE )

  #Construct the relevant A matrices and beta.shp
  A_list <- construct_Aobs_Ashp_betashp(yvec = yvec,
                                        mvec = mvec,
                                        ordering = ordering,
                                        inequalities_only = inequalities_only
                                        )

  A.shp <- A_list$A.shp
  A.obs <- A_list$A.obs
  beta.shp <- A_list$beta.shp

  ## Construct beta.obs and beta.obs_list ----

  #The following values are used by beta.obs
    #Note we need to do this here so that the yvalues remain constant across boostrap draws
  yvalues <- unique(yvec)
  mvalues <- unique(mvec)
  my_values <- purrr::cross_df(list(m=mvalues,y=yvalues)) %>%
    dplyr::arrange(m,y) %>%
    dplyr::select(y,m)

  #Compute beta.obs at the observed data
  beta.obs <- get_beta.obs_fn(yvec = yvec,
                           dvec = dvec,
                           mvec = mvec,
                           inequalities_only = inequalities_only,
                           yvalues = yvalues,
                           mvalues = mvalues,
                           my_values = my_values)

  #Bootstrap to get beta.obs_list
  beta.obs_list <- compute_bootstrap_draws_clustered(
    f =
      function(df,d,y,m,...){get_beta.obs_fn(
        yvec = df[[y]],
        dvec = df[[d]],
        mvec = df[[m]],
        inequalities_only = inequalities_only,
        yvalues = yvalues,
        mvalues = mvalues,
        my_values = my_values)},
    df = df,
    d = d,
    m = m,
    y = y,
    cluster = cluster,
    numdraws = B,
    return_df = F)

  ## Pass to the relevant moment inequality procedure ----

  if(method %in% c("FSST")){
    # Define target parameter
    A.tgt <- A_list$A.tgt

    # Run FSST
    lpm <- lpinfer::lpmodel(A.obs = A.obs,
                            A.shp = A.shp,
                            A.tgt = A.tgt,
                            beta.obs = beta.obs_list,
                            beta.shp = beta.shp)

    fsst_result <- lpinfer::fsst(df, lpmodel = lpm, beta.tgt = 0, R = B-1,
                                 weight.matrix = weight.matrix)

    return(list(result = fsst_result, reject = (fsst_result$pval < alpha)))
  }

  if(method == "CR"){
    stop("Cho and Russell not yet implemented")
  }

  if(method %in% c("ARP", "CS")){

    # Get variance matrix of the beta.obs boostraps
    sigma.obs <- stats::cov(base::Reduce(base::rbind,
                                         beta.obs_list))

    #Add the shape constraints as moments
    beta <- c(beta.obs, beta.shp)
    A <- rbind(A.obs, A.shp)

    #Update the covariance matrix to have zero in blocks for the shape constraints
    sigma <- rbind(cbind(sigma.obs, matrix(0,
                                           nrow = NROW(sigma.obs),
                                           ncol = NROW(A.shp) ) ),
                   matrix(0, nrow = NROW(A.shp), ncol = NCOL(sigma.obs) + NROW(A.shp) ))

    #Test if sigma is numerically not psd. If so, add small amount of noise
    min_eig <- base::min(base::eigen(sigma, only.values = TRUE)$values)
    if(min_eig < 0){
      sigma <- sigma + diag(10*abs(min_eig),
                            ncol = NROW(sigma),
                            nrow = NROW(sigma))}

    #Run the relevant test
    if(method == "ARP"){
      if(is.null(hybrid_kappa)){
      arp <- HonestDiD:::.lp_conditional_test_fn(theta = 0,
                                                 y_T = beta,
                                                 X_T = A,
                                                 sigma = sigma,
                                                 alpha = alpha,
                                                 hybrid_flag = "ARP")
      }else{
        lf_cv <- HonestDiD:::.compute_least_favorable_cv(X_T = A,
                                                         sigma = sigma,
                                                         hybrid_kappa = hybrid_kappa)

        hybrid_list <- list(hybrid_kappa = hybrid_kappa,
                            lf_cv = lf_cv)

        arp <- HonestDiD:::.lp_conditional_test_fn(theta = 0,
                                                   y_T = beta,
                                                   X_T = A,
                                                   sigma = sigma,
                                                   alpha = alpha,
                                                   hybrid_flag = "LF",
                                                   hybrid_list = hybrid_list)

      }

      return(arp)
    }

    if(method == "CS"){
      stop("Cox and Shi has not yet been implemented for the non-binary M case")
    }

  }


  stop("method must be one of ARP, CS, FSST, CR")

}


#Function for creating the matrices A.obs and A.shp and the vector beta.shp
construct_Aobs_Ashp_betashp <- function(yvec,
                                        mvec,
                                        ordering,
                                        inequalities_only = F){

  # Cardinality of the supports of Y and M
  d_y <- length(unique(yvec))
  K <- length(unique(mvec))

  # Make sure the ordering is defined for all support points
  if (!is.null(ordering) & !all(as.character(unique(mvec)) %in% names(ordering))) {
    stop("Variable ordering does not include all possible values of m!")
  }

  # Define "less than or equal to" function for the given partial ordering
  po_leq <- function(l, k) {
    ifelse(is.null(ordering), l <= k, l %in% ordering[[as.character(k)]])
  }

  # Getting ready to use lpinfer
  # x = (theta, delta, zeta, kappa, eta)
  # zeta, kappa: nuisance par to convert inequalities to equalities
  # eta: nuisance par to create target par (= theta_kk TV_kk)

  par_lengths <- c("theta" = K^2, "delta" = d_y * K, "zeta" = K,
                   "kappa" = d_y * K, "eta" = K)
  len_x <- sum(par_lengths)

  parameter_types <- c(rep("theta", par_lengths[["theta"]]),
                       rep("delta", par_lengths[["delta"]]),
                       rep("zeta", par_lengths[["zeta"]]),
                       rep("kappa", par_lengths[["kappa"]]),
                       rep("eta", par_lengths[["eta"]]))

  theta_indices <- which(parameter_types == "theta")
  delta_indices <- which(parameter_types == "delta")
  zeta_indices <- which(parameter_types == "zeta")
  kappa_indices <- which(parameter_types == "kappa")
  eta_indices <- which(parameter_types == "eta")



  # Set theta_lk = 0 for l > k using
  # Matrix that encodes whether l > k
  l_gt_k_mat <- matrix(NA, K, K)
  for (l in 1:K) {
    for (k in 1:K) {
      l_gt_k_mat[l, k] <- !po_leq(l, k)
    }
  }

  # Get hold of the indices
  l_gt_k_inds <- which(l_gt_k_mat)

  # Set shape constraints using A.shp x = beta.shp
  A.shp <- matrix(0, nrow = K, ncol = len_x)

  # Set bound on sum_{l!=k} theta_lk - sum_q delta_qk
  # Also includes eta_kk (= theta_kk TV_kk) as a "nuisance" parameter
  for (k in 1:K) {
    A.shp[k, ((k-1) * K + 1):(k * K)] <- 1
    A.shp[k, (k-1) * K + k] <- 0
    A.shp[k, par_lengths[1] + ((k-1) * d_y + 1):(k * d_y)] <- -1

    # Inequalities to equalities
    A.shp[k, sum(par_lengths[1:2]) + k] <- 1
    # TV_kk "nuisance" parameters
    A.shp[k, sum(par_lengths[1:4]) + k] <- 1
  }

  if(inequalities_only == T){
    #Remove both the extraneous thetas *and* kappa,zeta,eta
    A.shp <- A.shp[, -c(l_gt_k_inds, kappa_indices, eta_indices, zeta_indices)]

    #Add shape constraint that all parameters are >= 0 (this is not enforced by ARP)
    A.shp <- rbind(A.shp, diag(NCOL(A.shp)))

  }else{
    #Remove the extraneous thetas only
    A.shp <- A.shp[, -l_gt_k_inds]
  }

  beta.shp <- rep(0, NROW(A.shp))

  # Set remaining constraints using A.obs x = beta.obs

  # Define A.obs
  A.obs <- matrix(0,
                  nrow = K + K + (d_y * K),
                  ncol = len_x)

  for (k in 1:K) {
    # Match P(M = k | D = 0)
    A.obs[k, k + K * seq(0, K-1)] <- 1

    # Match P(M = k | D = 1)
    A.obs[K + k, ((k-1) * K + 1):(k * K)] <- 1

    ## # sum theta_{l<k} bound
    ## A.obs[2 * K + k, (k-1) * K + k] <- -1
    ## A.obs[2 * K + k, par_lengths[1] + ((k-1) * d_y + 1):(k * d_y)] <- -1

    # sup Delta(A) bound
    A.obs[2 * K + ((k-1) * d_y + 1):(k * d_y),
          par_lengths[1] + ((k-1) * d_y + 1):(k * d_y)] <- diag(d_y)

    # Inequalities to equalities
    ## A.obs[2 * K + k, sum(par_lengths[1:2]) + k] <- 1
    A.obs[cbind(2 * K + ((k-1) * d_y + 1):(k * d_y),
                sum(par_lengths[1:3]) + ((k-1) * d_y + 1):(k * d_y))] <- 1
  }

  if(inequalities_only){
  #Remove both the extraneous thetas *and* kappa,zeta,eta
  A.obs <- A.obs[, -c(l_gt_k_inds, kappa_indices, eta_indices, zeta_indices)]

  #The first 2K rows of A.obs are equality constraints, so we duplicate them with opposite signs to get equalities as inequalities
  A.obs <- rbind(A.obs[1:(2*K),],
                 -A.obs[1:(2*K),],
                 A.obs[(2*K+1):NROW(A.obs), ])
  }else{
    #Remove the extraneous thetas only
    A.obs <- A.obs[, -c(l_gt_k_inds)]
  }

  #Create A.tgt (only used for lpinfer functions)
  A.tgt <- numeric(len_x)
  A.tgt[sum(par_lengths[1:4]) + (1:K)] <- 1
  A.tgt <- A.tgt[-l_gt_k_inds]


  return(list(A.obs = A.obs,
              A.shp = A.shp,
              A.tgt = A.tgt,
              beta.shp = beta.shp))
}



# Constructing beta_obs
get_beta.obs_fn <- function(yvec, dvec, mvec, inequalities_only,
                            yvalues, mvalues, my_values) {
  #Get frequencies for all possible values of (y,m) | D=0
  p_ym_0_vec <- purrr::map_dbl(.x = 1:NROW(my_values),
                               .f = ~mean(yvec[dvec == 0] == my_values$y[.x]
                                          & mvec[dvec == 0] == my_values$m[.x]) )

  #Rearrange into a matrix where rows correspond to values of y and cols vals of m
  p_ym_d0 <-
    cbind(my_values, p_ym_0_vec) %>%
    reshape2::dcast(y ~ m, value.var = "p_ym_0_vec") %>%
    dplyr::select(-y) %>%
    as.matrix()


  #Get frequencies for all possible values of (y,m) | D=1
  p_ym_1_vec <- purrr::map_dbl(.x = 1:NROW(my_values),
                               .f = ~mean(yvec[dvec == 1] == my_values$y[.x]
                                          & mvec[dvec == 1] == my_values$m[.x]) )

  #Rearrange into a matrix where rows correspond to values of y and cols vals of m
  p_ym_d1 <-
    cbind(my_values, p_ym_1_vec) %>%
    reshape2::dcast(y ~ m, value.var = "p_ym_1_vec") %>%
    dplyr::select(-y) %>%
    as.matrix()

  # # Matrices (with dimension d_y x K) that store P(y,m | d) for d = 0,1
  # # d = 0
  # p_ym_d0 <- prop.table(table(factor(yvec)[dvec == 0], factor(mvec)[dvec == 0]))
  # # d = 1
  # p_ym_d1 <- prop.table(table(factor(yvec)[dvec == 1], factor(mvec)[dvec == 1]))

  # Matrices that store P(m | d) for d = 0,1
  # d = 0
  p_m_d0 <- colSums(p_ym_d0)
  # d = 1
  p_m_d1 <- colSums(p_ym_d1)

  if(inequalities_only){
    #Duplicate the first two sets of rows with opposite signs
    # to cast equalities as inequalities
    beta.obs <- c(p_m_d0, p_m_d1,
                  -p_m_d0, -p_m_d1,
                  p_ym_d1 - p_ym_d0)
  }else{
    beta.obs <- c(p_m_d0, p_m_d1,
                  p_ym_d1 - p_ym_d0)
  }

  return(beta.obs)
}
