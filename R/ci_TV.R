#' @title Construct confidence interval for total variation
#' @description This function calculates a confidence interval for the total
#'   variation for a given always taker group. The inference is via inverting
#'   the test of Fang et al. (2023; FSST)
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable, which is assumed to take a discrete
#'   support
#' @param at_group Indicates for which always-taker group the confidence
#'   interval is calculated for. It can take any value in the support of M.
#' @param ordering A list with length equal to the cardinality of the support of
#'   the mediator variable. The name of each element corresponds to a point in
#'   the support, and each element is a vector that collects all m values that
#'   are less than or equal to this point. If ordering = NULL, the standard
#'   ordering is used.
#' @param B Bootstrap size, default is zero
#' @param alpha Significance level
#' @param weight.matrix Weight matrix used to implement FSST. Possible options
#'   are "diag", "avar", "identity." Defaults is "diag" as in FSST.
#' @param bisec If bisec=TRUE, a bisection method is used rather than a grid
#'   search
#' @param verbose Prints test results for each null value if verbose = TRUE
#' @param eps Tolerance level for the bisection method (in terms of the
#'   p-value), default is 10\% of the significance level alpha
#' @export
ci_TV <- function(df,
                  d,
                  m,
                  y,
                  at_group,
                  ordering = NULL,
                  B = 500,
                  alpha = .05,
                  weight.matrix = "diag",
                  bisec = FALSE,
                  verbose = FALSE,
                  eps = .1 * alpha){

  yvec <- df[[y]]
  dvec <- df[[d]]
  mvec <- df[[m]]

  # Cardinality of the supports of Y and M
  d_y <- length(unique(yvec))
  K <- length(unique(mvec))

  # Check at_group actually exists
  if (!(at_group %in% unique(mvec))) {
    stop("at_group must be a support point of M!")
  }
  
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

  # Set theta_lk = 0 for l > k using
  # Matrix that encodes whether l > k
  l_gt_k_mat <- matrix(NA, K, K)
  for (l in 1:K) {
    for (k in 1:K) {
      l_gt_k_mat[l, k] <- !po_leq(l, k)
    }
  }

  # Get hold of the indices; drop these parameters later
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
    A.shp[k, sum(par_lengths[1:2]) + k] <- -1
    # TV_kk "nuisance" parameters
    A.shp[k, sum(par_lengths[1:4]) + k] <- 1
  }
  A.shp <- A.shp[, -l_gt_k_inds]
  beta.shp <- rep(0, K)


  # Set remaining constraints using A.obs x = beta.obs
  # Bootstrap sample indices
  boot_mat <- matrix(sample(1:nrow(df), B * nrow(df), replace = T), nrow = B)  
  # Get all beta_obs to pass to lpinfer
  beta.obs_list <- lapply(1:B,
                          function(b) get_beta.obs(factor(yvec)[boot_mat[b,]],
                                                   dvec[boot_mat[b,]],
                                                   factor(mvec)[boot_mat[b,]])
                          )

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
                sum(par_lengths[1:3]) + ((k-1) * d_y + 1):(k * d_y))] <- -1
  }

  A.obs <- A.obs[,-l_gt_k_inds]

  # Function to test H_0: TV_kk = t
  test_TV <- function(t) {
    # Define target parameter
    A.tgt <- numeric(len_x)
    A.tgt[sum(par_lengths[1:4]) + at_group] <- 1
    A.tgt[(at_group - 1) * K + at_group] <-  -t
    
    A.tgt <- A.tgt[-l_gt_k_inds]

    # Run FSST
    lpm <- lpinfer::lpmodel(A.obs = A.obs,
                            A.shp = A.shp,
                            A.tgt = A.tgt,
                            beta.obs = beta.obs_list,
                            beta.shp = beta.shp)

    result <- lpinfer::fsst(df,
                            lpmodel = lpm,
                            beta.tgt = 0,
                            R = B-1,
                            weight.matrix = weight.matrix)

    return(result)
  }
  
  if (!bisec) {
    # Take a grid on null values
    # Change later to a smarter search method
    tgrid <- seq(0, 1, .02)

    # vector to collect all p-values
    p_vals <- numeric(length(tgrid))

    # Test!
    for (t in tgrid) {
      print(paste0("Working on null value = ", t))
      j
      result <- test_TV(t)
      
      if (verbose) {
        print(c(t, result$pval[2]))
      }    
      
      p_vals[which(t == tgrid)] <- result$pval[2]
    }
  } else {


    if ((t0_pval <- test_TV(0)$pval[2]) > alpha) {
      
      return(cbind(0, t0_pval))
      
    } else {
      
      tgrid <- NULL
      t <- .5
      p_vals <- NULL
      pval <- 0

      # Endpoints for bisection method
      a <- 0
      b <- 1
      
      while(abs(pval - .05) > eps) {

        print(paste0("Working on null value = ", t))
        
        tgrid <- c(tgrid, t)
        result <- test_TV(t)

        p_val <- result$pval[2]
        p_vals <- c(p_vals, pval)

        if (verbose) {
          print(c(t, p_val))
        }
        # Update null value
        if (p_val > .05) {

          b <- t
          t <- (t + a) / 2
          
        } else {
          
          a <- t
          t <- (t + b) / 2
          
        }
      }   
    }
  }
  
  return(cbind(tgrid, p_vals))  
}
