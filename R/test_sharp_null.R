#' @title Hypothesis test for the sharp null of full mediation
#' @description This function tests the sharp null that Y(1,m) = Y(0,m). The
#'   outcome and mediator are both assumed to take finitely many different
#'   values. Several options are available for infernece.
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable(s)
#' @param y Name of the outcome variable, which is assumed to take a discrete
#'   support
#' @param method Method to use. One of "ARP, "CS", "FSST", "CR"
#' @param ordering A list with length equal to the cardinality of the support of the mediator variable. The name of each element corresponds to a point in the support, and each element is a vector that collects all m values that are less than or equal to this point. If ordering = NULL, the standard ordering is used. If length(m) > 1, then the default is the elementwise ordering.
#' @param B Bootstrap size, default is 500
#' @param cluster Cluster for bootstrap
#' @param weight.matrix Weight matrix used to implement FSST. Possible options are "diag", "avar", "identity." Defaults is "diag" as in FSST.
#' @param hybrid_kappa The first-stage size value of the ARP hybrid test. If NULL, the ARP conditional (non-hybrid) test is used. Default is alpha/10
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
#' @param alpha Significance level. Default is 0.05
#' @param rearrange Logical variable indicating whether to rearrange the conditional probabilities to obey monotonicity. De
#' @param eps_bar Cho and Russell (2023) truncation parameter
#' @param enumerate Enumerate vertices for Cox and Shi (2023) implementataion?
#' @param fix_n1 Whether the number of treated units (or clusters) should be fixed in the bootstrap
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
                            alpha = 0.05,
                            rearrange = FALSE,
                            eps_bar = 1e-03,
                            enumerate = FALSE,
                            fix_n1 = TRUE){

  ## Process the inputted df ----

  # Remove missing values

  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)

  # Sample size
  n <- nrow(df)

  # If M is multivariate convert to univariate and use elementwise ordering
  if (length(m) > 1) {

    m_supp <- unique(df[m])
    ordering <- vector("list", length = nrow(m_supp))
    names(ordering) <- 1:nrow(m_supp)

    # Create ordering and replace multivariate M by a univariate M
    uni_m <- numeric(n)
    for (k in 1:nrow(m_supp)) {
      m_k <- m_supp[k,]
      ordering[[as.character(k)]] <-
        which(apply(m_supp, 1, function(m_supp_elem) all(m_supp_elem <= m_k)))
      uni_m[colSums(as.vector(m_k) == t(df[m])) == 2] <- k
    }

    df$m <- uni_m
    m <- "m"
  }

  #Discretize y if needed
  if(!is.null(num_Ybins)){
    df[[y]] <- discretize_y(yvec = df[[y]], numBins = num_Ybins)
  }

  yvec <- df[[y]]
  dvec <- df[[d]]
  mvec <- df[[m]]

  ## Construct the A matrices and beta.shp ----

  #Specify whether method requires us to input only inequalities
  inequalities_only <- ifelse(method %in% c("ARP","CS", "CR"),
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

  # Current version allows only rearrange = TRUE when method = CR
  if (method == "CR") {
    regarrange <- TRUE
  }

  # Compute beta.obs at the observed data
  beta.obs <- get_beta.obs_fn(yvec = yvec,
                              dvec = dvec,
                              mvec = mvec,
                              inequalities_only = inequalities_only,
                              yvalues = yvalues,
                              mvalues = mvalues,
                              my_values = my_values,
                              rearrange = rearrange)

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
        my_values = my_values,
        rearrange = rearrange)},
    df = df,
    d = d,
    m = m,
    y = y,
    cluster = cluster,
    numdraws = B,
    return_df = F,
    fix_n1 = fix_n1)

  ## Pass to the relevant moment inequality procedure ----

  if(method %in% c("FSST")){
    # Define target parameter
    A.tgt <- A_list$A.tgt

    # Run FSST
    lpm <- lpinfer::lpmodel(A.obs = A.obs,
                            A.shp = A.shp,
                            A.tgt = A.tgt,
                            beta.obs = c(list(beta.obs),
                                         beta.obs_list),
                            beta.shp = beta.shp)

    fsst_result <- lpinfer::fsst(df, lpmodel = lpm, beta.tgt = 0, R = B,
                                 weight.matrix = weight.matrix)

    return(list(result = fsst_result, reject = (fsst_result$pval[1, 2] < alpha)))
  }

  if (method == "CR") {
    # Define target parameter
    A.tgt <- A_list$A.tgt
    len_x <- length(A.tgt)

    params <- list(OutputFlag=0)


    # Define gurobi model
    model <- list()

    A <- rbind(A.shp, A.obs)
    rhs <- c(beta.shp, beta.obs)
    lb <- rep(0, len_x)
    ub <- rep(1, len_x)

    # Combine lower and upper bound into matrix
    model$A <- A
    model$obj <- A.tgt
    model$rhs <- rhs
    model$lb  <- lb
    model$ub  <- ub
    model$sense <- rep('>', length(rhs))

    # Optimize for "l"
    model$modelsense <- 'min'
    min.result <- gurobi::gurobi(model, params)

    # Optimize for "u"
    model$modelsense <- 'max'
    max.result <- gurobi::gurobi(model, params)

    # Draw perturbation
    xi_obj <- runif(len_x) * eps_bar
    xi_rhs <- runif(length(rhs)) * eps_bar
    xi_lb <- runif(len_x) * eps_bar
    xi_ub <- runif(len_x) * eps_bar

    ############################################################################
    # Compute LB-,LB+,UB-,UB+
    ############################################################################

    model$rhs <- rhs - xi_rhs
    model$lb <- lb - xi_lb
    model$ub <- ub + xi_ub

    ############################################################################
    # Compute LB-,UB-
    model$obj<- A.tgt - xi_obj

    # Optimize LB-
    model$modelsense <- 'min'
    min.result.m <- gurobi::gurobi(model, params)

    # Record lower bound
    lbminus <- min.result.m$objval

    # Optimize UB-
    model$modelsense <- 'max'
    max.result.m <- gurobi::gurobi(model, params)

    # Record lower bound
    ubminus <- max.result.m$objval

    #######################################################
    # Compute LB+,UB+
    model$obj <- A.tgt + xi_obj

    # Optimize
    model$modelsense <- 'min'
    min.result.p <- gurobi::gurobi(model, params)

    # Record lower bound
    lbplus <- min.result.p$objval

    # Optimize UB-
    model$modelsense <- 'max'
    max.result.p <- gurobi::gurobi(model, params)

    # Rrecord upper bound
    ubplus <- max.result.p$objval

    ############################################################################
    # Begin bootstrap procedure
    boot_lbminus <- rep(NA, B)
    boot_lbplus <- rep(NA, B)
    boot_ubminus <- rep(NA, B)
    boot_ubplus <- rep(NA, B)

    for (b in 1:B) {

      # Get beta.obs for bootstrap draw
      beta.obs_b <- beta.obs_list[[b]]

      # Update rhs
      rhs_b <- c(beta.shp, beta.obs_b)

      # Update model
      model$rhs <- rhs_b - xi_rhs

      ############################################################################
      # Compute LB-,UB-
      model$obj<- A.tgt - xi_obj

      # Optimize LB-
      model$modelsense <- 'min'
      bmin.result.m <- gurobi::gurobi(model, params)

      # Record lower bound
      boot_lbminus[b] <- bmin.result.m$objval

      # Optimize UB-
      model$modelsense <- 'max'
      bmax.result.m <- gurobi::gurobi(model, params)

      # Record upper bound
      boot_ubminus[b] <- bmax.result.m$objval

      ###################################################
      # Objective function for LB+ and UB+
      model$obj<- A.tgt + xi_obj

      # Optimize LB+
      model$modelsense <- 'min'
      bmin.result.p <- gurobi::gurobi(model, params)

      # Record lower bound
      boot_lbplus[b] <- bmin.result.p$objval

      # Optimize UB+
      model$modelsense <- 'max'
      bmax.result.p <- gurobi::gurobi(model, params)

      # Record upper bound
      boot_ubplus[b] <- bmax.result.p$objval
    }

    # Compute indicator Dn
    bn <- 1/sqrt(log(n))
    Delta <- max(ubplus, ubminus) - min(lbminus, lbplus)
    Dn <- (Delta > bn) + 0

    # Calculate kappa
    kappa <- (1 - alpha) * Dn + (1 - alpha/2) * (1 - Dn)

    # Bootstrap quantities
    lbminus_q <- sqrt(n) * (boot_lbminus - lbminus)
    lbplus_q <- sqrt(n) * (boot_lbplus - lbplus)
    ubminus_q <- -sqrt(n) * (boot_ubminus - ubminus)
    ubplus_q <- -sqrt(n) * (boot_ubplus - ubplus)

    # Select quantile according to kappa
    psi_k_lb_minus <- quantile(lbminus_q, kappa)
    psi_k_lb_plus <- quantile(lbplus_q, kappa)
    psi_k_ub_minus <- quantile(ubminus_q, kappa)
    psi_k_ub_plus <- quantile(ubplus_q, kappa)

    #compute confidence set for alpha=0.05
    CSlb <- min(lbminus, lbplus) - (1/sqrt(n)) * max(psi_k_lb_minus, psi_k_lb_plus)
    CSub <- max(ubminus, ubplus) + (1/sqrt(n)) * max(psi_k_ub_minus, psi_k_ub_plus)

    return(list(CI = c(CSlb, CSub), reject = (0 < CSlb)))
  }

  if(method %in% c("ARP", "CS")){

    # Get variance matrix of the beta.obs boostraps
    sigma.obs <- stats::cov(base::Reduce(base::rbind,
                                         beta.obs_list))

    #Add the shape constraints as moments
    beta <- c(beta.obs, beta.shp)
    A <- rbind(A.obs, A.shp)


    #Run the relevant test
    if (method == "ARP" ) {

      #Update the covariance matrix to have zero in blocks for the shape constraints
      sigma <- rbind(cbind(sigma.obs, matrix(0,
                                             nrow = NROW(sigma.obs),
                                             ncol = NROW(A.shp) ) ),
                     matrix(0, nrow = NROW(A.shp), ncol = NCOL(sigma.obs) + NROW(A.shp) ))

      #Test if sigma is numerically not psd. If so, add small amount of noise
      min_eig <- base::min(base::eigen(sigma, only.values = TRUE)$values)
      if (min_eig < 0) {
        sigma <- sigma + diag(10 * abs(min_eig),
                              ncol = NROW(sigma),
                              nrow = NROW(sigma))
      }

      if (is.null(hybrid_kappa)) {
      arp <- HonestDiD:::.lp_conditional_test_fn(theta = 0,
                                                 y_T = beta,
                                                 X_T = A,
                                                 sigma = sigma,
                                                 alpha = alpha,
                                                 hybrid_flag = "ARP")
      } else {
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

      sigma <- sigma.obs

      ## Convert in C * delta <= m for to make use of the replication package
      C_Z <- A
      B_Z <- rbind(diag(length(beta.obs)),
                   matrix(0, nrow = length(beta.shp), ncol = length(beta.obs)))
      d_Z <- 0

      if (min(base::eigen(sigma, only.values = T)$values) < 10^-6){
        #If sigma is not full-rank we extract the full rank component

        eigendecomp <- eigen(sigma)
        eigenvals <- eigendecomp$values
        eigenvecs <- eigendecomp$vectors

        tol <- 10^-6
        positive_indices <- which(eigenvals > tol)

        #Create a matrix that selects rows with positive indices
        rem_red <- diag(length(eigenvals))
        rem_red <- as.matrix(rem_red[positive_indices,])

        #Deal with annoying fact that R converts to column vector if length 1
        if (length(positive_indices) == 1) { rem_red <- t(rem_red) }

        rem_red <- eigenvecs %*% t(rem_red)

        beta.obs_red <- t(rem_red) %*% beta.obs
        V_red <- diag( x = eigenvals[positive_indices],
                        nrow = length(eigenvals[positive_indices]))


        # Adjust the linear programming parameters accordingly
        # We have B_Z_red %*% mu - C_Z delta <= d_Z, following CS notation
        B_Z_red <- B_Z %*% rem_red
        d_Z <- B_Z_red %*% beta.obs_red - B_Z %*% beta.obs
        sigma <- V_red
        sigmaInv <- solve(sigma)
      } else {
        beta.obs_red <- beta.obs
        B_Z_red <- B_Z
        sigmaInv <- solve(sigma)
      }

      d_nuis <- ncol(A)  # number of nuisance parameters
      d_ineq <- nrow(A)  # number of inequalities

      # Perform CC Test

      #Convert beta.obs_red to a column vector if not already
      beta.obs_red <- as.matrix(beta.obs_red, ncol = length(beta.obs_red))

      # Finding T (Test Statistic)

      # Amat x <= d_Z
      Amat <- rbind(cbind(B_Z_red, - C_Z),
                    cbind(matrix(0, nrow = d_nuis, ncol = ncol(B_Z_red)), diag(d_nuis)))
      u <- c(d_Z, rep(1, d_nuis))
      l <- c(rep(-Inf, length(d_Z)), rep(0, d_nuis))
      Dmat <- matrix(0, nrow = (nrow(sigma) + d_nuis), ncol = (nrow(sigma) + d_nuis))
      Dmat[1:nrow(sigma), 1:nrow(sigma)] <- 2 * sigmaInv
      dvec <- matrix(c(2 * sigmaInv %*% beta.obs_red, numeric(d_nuis)), ncol = 1)

      qp <- osqp::solve_osqp(P = Dmat, q = -dvec, A = Amat, l = l, u = u, list(verbose = FALSE))
      ## qp <- quadprog::solve.QP(Dmat = Dmat,
      ##                          dvec = dvec,
      ##                          Amat = Amat,
      ##                          bvec = - d_Z)



      ## beta.obs_red_star <- qp$solution[1:nrow(sigma)]
      beta.obs_red_star <- qp$x[1:nrow(sigma)]
      T_CC <- t(beta.obs_red - beta.obs_red_star) %*%
        sigmaInv %*%
        (beta.obs_red - beta.obs_red_star)


      # Find the Degree of Freedom:
      tol <- 1e-10

      if (enumerate) {
        A_vert <- rbind(-diag(d_ineq),
                        t(C_Z),
                        -t(C_Z),
                        rep(1, d_ineq),
                        -rep(1, d_ineq))
        b_vert <- c(numeric(d_ineq), 2 * numeric(d_nuis), 1, -1)
        H <- vertexenum::enumerate.vertices(A = A_vert, b = b_vert)
        A_Z <- H %*% B_Z_red
        b_Z <- H %*% d_Z

        # Binding constraints

        dof_n <- sum(abs(A_Z %*% beta.obs_red_star - b_Z) > tol)
      } else {
        # Define A_ineq0 and b_ineq0
        A_ineq0 <- diag(d_ineq)
        ## b_ineq0 <- matrix(0, nrow = d_ineq, ncol = 1)

        # Define A_eq0, A_eq, and b_eq
        # mstar = - d_Z + B_Z_red %*% beta.obs_red*
        mstar <- - d_Z + B_Z_red %*% beta.obs_red_star
        A_eq0 <- rbind(t(C_Z), t(mstar), rep(1, d_ineq))
        A_eq <- rbind(t(C_Z), rep(1, d_ineq))
        b_eq <- c(rep(0, d_nuis), 1)

        # Perform linear programming to find Vmu_min
        ## lp_options <- lpSolve::lp.control(display = "silent", simplex = "dual")
        ## result <- lpSolve::lp("min", mstar, rbind(A_ineq0, A_eq),
        ##                       c(rep(">=", d_ineq), rep("=", nrow(A_eq))),
        ##                       c(b_ineq0, b_eq))
        ## Vmu_min <- result$objval
        ## flag <- result$status

        model <- list()
        model$A <- A_eq
        model$obj <- mstar
        model$modelsense <- "min"
        model$rhs <- b_eq
        model$sense <- rep("=", nrow(A_eq))
        model$lb <- rep(0, d_ineq)

        result <- gurobi::gurobi(model)

        Vmu_min <- result$objval
        flag <- result$status

        # Check conditions and calculate dof_n
        if (flag != "OPTIMAL") {
          dof_n <- 0
          ## } else if (flag == 2) { ## for lpSolve
        } else if (Vmu_min >= 0.00005) {
          dof_n <- 0
        } else {
          # Initialize nb_ineq_min
          nb_ineq_min <- rep(NA, d_ineq)
          counter <- 1

          # Iterate through each j
          for (bj in 1:d_ineq) {
            # Find the largest b_j allowed

            ## result_j <- lpSolve::lp("min", A_ineq0[bj,],
            ##                         rbind(A_eq0, diag(d_ineq), diag(d_ineq)),
            ##                         c(rep("=", nrow(A_eq0)), rep(">=", d_ineq), rep("<=", d_ineq)),
            ##                         c(rep(0, d_nuis), Vmu_min, 1, b_ineq0, b_ineq0 + 1))

            model <- list()
            model$A <- A_eq0
            model$obj <- -A_ineq0[bj, ]
            model$modelsense <- "min"
            model$rhs <- c(rep(0, d_nuis), Vmu_min, 1)
            model$sense <- rep("=", nrow(A_eq0))
            model$lb <-  0
            model$ub <- 1

            result_j <- gurobi::gurobi(model)
            # Update nb_ineq_min
            ## if (result_j$status == 2) {# for lpSolve
            if (result_j$status != "OPTIMAL") {
              nb_ineq_min[bj] <- - 1
              counter <- counter + 1
            } else {
              nb_ineq_min[bj] <- result_j$objval
            }
          }

          # Collect rows of A_ineq corresponding to implicit equalities
          A_impeq <- A_ineq0[(0 - nb_ineq_min) < tol, ]

          # Combine A_impeq and A_eq0
          A_full_eq <- rbind(A_impeq, A_eq0)

          # Calculate rank of A_full_eq
          rkA <- qr(A_full_eq)$rank

          # Calculate dof_n
          dof_n <- d_ineq - rkA + 1
        }

      }
      # Enumerating the vertices using (outdated) package vertexenum
      # Finds vertices for the polytope Ax <= b

      cv_CC <- qchisq(1 - alpha, df = dof_n)
      return(list(T_CC = T_CC, cv_CC = cv_CC, df = dof_n, reject = (T_CC > cv_CC)))
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
    A.shp[k, sum(par_lengths[1:2]) + k] <- -1
    # TV_kk "nuisance" parameters
    A.shp[k, sum(par_lengths[1:4]) + k] <- 1
  }

  if (inequalities_only == T) {
    #Remove both the extraneous thetas *and* kappa,zeta,eta
    A.shp <- A.shp[, -c(l_gt_k_inds, kappa_indices, eta_indices, zeta_indices)]

    #Add shape constraint that all parameters are >= 0 (this is not enforced by ARP)
    A.shp <- rbind(A.shp, diag(NCOL(A.shp)))

  } else {
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
                sum(par_lengths[1:3]) + ((k-1) * d_y + 1):(k * d_y))] <- -1
  }

  if (inequalities_only) {
  #Remove both the extraneous thetas *and* kappa,zeta,eta
  A.obs <- A.obs[, -c(l_gt_k_inds, kappa_indices, eta_indices, zeta_indices)]

  #The first 2K rows of A.obs are equality constraints, so we duplicate them with opposite signs to get equalities as inequalities
  A.obs <- rbind(A.obs[1:(2 * K),],
                 -A.obs[1:(2 * K),],
                 A.obs[(2 * K + 1):NROW(A.obs), ])
  } else {
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
              beta.shp = beta.shp,
              par_length = length(A.tgt)))
}



# Constructing beta_obs
get_beta.obs_fn <- function(yvec, dvec, mvec, inequalities_only,
                            yvalues, mvalues, my_values, rearrange = FALSE) {
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

  if (rearrange) {

    cdf_m_d0 <- pmax(cumsum(p_m_d0), cumsum(p_m_d1))
    cdf_m_d1 <- pmin(cumsum(p_m_d0), cumsum(p_m_d1))

    p_m_d0 <- diff(c(0, cdf_m_d0))
    p_m_d1 <- diff(c(0, cdf_m_d1))

  }

  if (inequalities_only) {
    #Duplicate the first two sets of rows with opposite signs
    # to cast equalities as inequalities
    beta.obs <- c(p_m_d0, p_m_d1,
                  -p_m_d0, -p_m_d1,
                  p_ym_d1 - p_ym_d0)
  } else {
    beta.obs <- c(p_m_d0, p_m_d1,
                  p_ym_d1 - p_ym_d0)
  }

  return(beta.obs)
}
