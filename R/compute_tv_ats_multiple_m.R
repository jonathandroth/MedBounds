#' @title Bounds on fraction of always-takers affected with a multi-dimensional m
#'@description Computes  a lower bound on the TV distance btwn Y(1,m) and Y(0,m) for always takers
#' who have M(1)=M(0)=m
#' This is a lower bound on the fraction of these ATs for whom there is a direct effect of D on Y
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Vector of the mediator variable
#' @param y Name of the outcome variable
#' @param w (Optional) Name of weighting variable. If null, equal weights are used
#' @param at_group (Optional) Value of m specifying which always-takers to compute lower bounds of TV for.
#' If at_group is specified, then we compute a lower bound on TV between Y(1,at_group) and Y(0,at_group) for
#' ATs who have M(1)=M(0)=at_group. If at_group is null (the default), we compute a lower bound on
#' the weighted average of TV across all always-takers, with weights proportional to shares in population
#' @export

compute_tv_ats_multiple_m <- function(df,
                                      d,
                                      m,
                                      y,
                                      at_group = NULL,
                                      w = NULL){

  yvec <- df[[y]]
  dvec <- df[[d]]
  mdf <- df[,m]

  if(is.null(w)){
    wvec <- rep(1, NROW(df))
  }else{
    wvec <- df[[w]]
  }

  max_p_diffs_list <- compute_max_p_difference(dvec = dvec,
                                               mdf = mdf,
                                               yvec = yvec,
                                               wvec=wvec)
  max_p_diffs <- max_p_diffs_list$max_p_diffs
  mvalues <- max_p_diffs_list$mvalues


    #Now, we create a linear program to find compliance types that minimize the TV distance

    #Create a list of all possible compliance types. This is all possible pairwise combinations of unique rows of mdf
     #First, create a data frame containing all combinations of the row numbers of mvalues
    mvalues_df <- base::expand.grid(1:NROW(mvalues), 1:NROW(mvalues))
    m0_types <- mvalues[mvalues_df[,1],]
    m1_types <- mvalues[mvalues_df[,2],]

    #Determine which rows of m0_types are weakly less than the corresponding row in m1_types
    monotonic_types <- base::sapply(1:NROW(m0_types),
                                    FUN = function(i){
                                      base::all(m0_types[i,] <= m1_types[i,])})

    #Restrict to the types satisfying monotonicity
    m0_types <- m0_types[monotonic_types,]
    m1_types <- m1_types[monotonic_types,]


    #Compute the marginal distribution of M among D=1
    p_m_1_fn <- function(mvalue){
      mdf1 <- mdf[dvec ==1,]
      M_equals_mvalue <- sapply(1:NROW(mdf1), function(i){all(mdf1[i,]==mvalue)} )
      return(stats::weighted.mean( x = M_equals_mvalue,
                                   weights = wvec[dvec == 1]/sum(wvec[dvec == 1])))
    }

    p_m_1 <- base::apply(mvalues,1, p_m_1_fn)

    #Compute the marginal distribution of M among D=0
    p_m_0_fn <- function(mvalue){
      mdf0 <- mdf[dvec ==0,]
      M_equals_mvalue <- sapply(1:NROW(mdf0), function(i){all(mdf0[i,]==mvalue)} )
      return(stats::weighted.mean( x = M_equals_mvalue,
                                   weights = wvec[dvec == 0]/sum(wvec[dvec == 0])))
    }

    p_m_0 <- base::apply(mvalues,1, p_m_0_fn)

    ## We now estimate the TV_kk
    ## We consider an optimization where the first part of the optimization vector corresponds to the shares of complier types
    ## and the final vector is the violation of the moments (corresponding to our lb on theta_kk TV_kk)

    ##Create a matrix corresponding to the constraint that the sum of the m1_types must equal the marginal distribution of M among D=1
    #Create a matrix where each row i is a NROW(m1_types) length vector where the jth element indicates if the jth row of m1_types equals the ith row of mvalues
    m1_marginals_constraints_matrix <-
    Reduce(rbind,
           base::lapply(X = 1:NROW(mvalues),
                 FUN = function(i){base::apply(m1_types, 1, function(x){ base::all(x == mvalues[i,]) })}))

    #Add a matrix of zeros to m1_marignals_constraints_matrix with length NROW(mvalues)
    m1_marginals_constraints_matrix <- base::cbind(m1_marginals_constraints_matrix,
                                                   matrix(0,
                                                          nrow = NROW(m1_marginals_constraints_matrix),
                                                          ncol = NROW(mvalues)))

    ##Create a matrix corresponding to the constraint that the sum of the m1_types must equal the marginal distribution of M among D=0
    #Create a matrix where each row i is a NROW(m0_types) length vector where the jth element indicates if the jth row of m0_types equals the ith row of mvalues
    m0_marginals_constraints_matrix <-
      Reduce(rbind,
             base::lapply(X = 1:NROW(mvalues),
                  FUN = function(i){base::apply(m0_types, 1, function(x){ base::all(x == mvalues[i,]) })}))

    #Add a matrix of zeros to m0_marignals_constraints_matrix
    m0_marginals_constraints_matrix <- base::cbind(m0_marginals_constraints_matrix,
                                                   matrix(0,
                                                          nrow = NROW(m0_marginals_constraints_matrix),
                                                          ncol = NROW(mvalues)))


    ##Create a constraint corresponding to the constraint that max_p_diffs must be greater than or equal to the sum of all complier types that end up at a given mvalue
    #Create a matrix where each row i is a NROW(m1_types) length vector where the jth element indicates if the jth row of m1_types equals the ith row of mvalues and the jth row of m1_types is not exactly equal to the jth row of m0_types
    complier_constraints_matrix <-
      Reduce(rbind,
        base::lapply(X = 1:NROW(mvalues),
                     FUN = function(i){base::sapply(1:NROW(m1_types), function(s){ base::all(m1_types[s,] == mvalues[i,]) & !base::all(m0_types[s,] == m0_types[i,]) })}))

    #Add an identity matrix to complier_constraints_matrix
      # Thus the end vector corresponds to violation of kth moment, i.e. theta_kk TV_k
    complier_constraints_matrix <- base::cbind(complier_constraints_matrix,
                                               diag(NROW(mvalues)))

    #Combine the constraint matrices
    constraints_matrix <- base::rbind(m1_marginals_constraints_matrix, m0_marginals_constraints_matrix, complier_constraints_matrix)

    #Combine the constants associated with the matrices
    d <- c(p_m_1, p_m_0, max_p_diffs)

    #Specify the direction of the equalities/inequalities
    dir <- c(rep("==", 2*NROW(m1_marginals_constraints_matrix)), rep(">=", NROW(max_p_diffs)))

    #Specify the bounds on the optimization vector:
    # The elements in the first part of the optimizaiton vector are probabilities, thus btwn 0 and 1
    # The bounds theta_kk TV_k must be >= 0
    # Put this in a list with elements lower and upper
    bounds <- list(lower = list(ind = 1:NCOL(constraints_matrix),
                                val = rep(0, NCOL(constraints_matrix))),
                   upper = list(ind = 1:NCOL(constraints_matrix),
                                val = c(rep(1, NCOL(constraints_matrix)-NROW(mvalues)),
                                        rep(Inf, NROW(mvalues))) )
                    )

    ##Specify the objective

    ##Create a function that returns the objective vectors for a single at_group
    ##Our objective is theta_kk TV_k / theta_kk
    compute_obj_vectors <- function(at_group){

      row_equals <- function(x,y){ base::all(x == y)}

      #Find the index i such that the ith row of m0_types and m1_types are equal to at_group
      at_group_index <- which(purrr::map_lgl(.x = 1:NROW(m0_types),
                                             .f = ~row_equals(m0_types[.x,], at_group) &
                                               row_equals(m1_types[.x,], at_group)))

      #Thus denominator is a basis vector with one in pos corresponding to thetakk
      obj_denominator <- rep(0, NCOL(constraints_matrix))
      obj_denominator[at_group_index] <- 1

      #Find the index i such that the ith row of mvalues equals at_group
      at_group_index <- which(purrr::map_lgl(.x = 1:NROW(mvalues),
                                             .f = ~row_equals(mvalues[.x,],at_group)))

      #Numerator is a basis vector with one in pos corresponding to thetakk TVk
      obj_numerator <- rep(0, NCOL(constraints_matrix))
      obj_numerator[NROW(m1_types) + at_group_index] <- 1

      return(list(obj_numerator = obj_numerator,
                  obj_denominator = obj_denominator))
    }

    if(!is.null(at_group)){
      obj_numerator <- compute_obj_vectors(at_group = at_group)$obj_numerator
      obj_denominator <- compute_obj_vectors(at_group = at_group)$obj_denominator
    }else{
      #Otherwise, we're interested in the average, i.e. \sum_k \theta_kk TV_k / \sum_k \theta_kk
      #Hence, we compuute the objective vector for each at_group and take the average
      obj_numerator <- Reduce("+", base::lapply(X = 1:NROW(mvalues),
                                                FUN = function(i){compute_obj_vectors(at_group = mvalues[i,])$obj_numerator}))
      obj_denominator <- Reduce("+", base::lapply(X = 1:NROW(mvalues),
                                                 FUN = function(i){compute_obj_vectors(at_group = mvalues[i,])$obj_denominator}))
    }
    #Use Rglpk to solve the fraction-linear program to minimize the weighted average of TV
      #We capture warnings if the denom can be zero
    quiet_fractional_LP <- purrr::quietly(MedBounds::Rglpk_solve_fractional_LP)
    max_violation <-
      quiet_fractional_LP(
        obj_numerator = obj_numerator,
        obj_denominator = obj_denominator,
        mat = constraints_matrix,
        dir = dir,
        rhs = d,
        bounds = bounds,
        max = FALSE)

    if(!is.na(max_violation$warnings[1]) & max_violation$warnings[1] == "The minimum value of the denominator is not positive. Returning NaN."){
      warning("The lower bound for the fraction of always-takers for the chosen at_group is zero. Returning NaN")
      return(NaN)
    }
    #Return the maximal violation
    return(max_violation$result$optimum)
}


#' @title Wrapper for Rglpk_solve_lp that implements linear fractional programming
#' @param obj_numerator The objective vector for the numerator of the linear fractional program
#' @param obj_deonminator The objective vector for the denominator of the linear fractional program. Note that the denominator must be strictly positive over all feasible values (if not, we return a warning and NaN for the value)
#' @param mat The constraint matrix
#' @param dir The direction of the constraints
#' @param rhs The right hand side of the constraints
#' @param bounds The bounds on the optimization vector
#' @param max (Optional) Whether to maximize or minimize the objective
#' @param constant_numerator (Optional) A constant to add to the numerator of the linear fractional program
#' @param constant_denominator (Optional) A constant to add to the denominator of the linear fractional program
#' @param ... (Optional) Additional arguments to be passed to Rglpk_solve_lp
#' @return The optimal value of the linear fractional program

Rglpk_solve_fractional_LP <- function(obj_numerator,
                                      obj_denominator,
                                      mat,
                                      dir,
                                      rhs,
                                      bounds,
                                      max = FALSE,
                                      constant_numerator = 0,
                                      constant_denominator = 0,
                                      ...){

  #First, check that the minimum value of the denominator is positive
  denom_lp <- Rglpk::Rglpk_solve_LP(obj = obj_denominator,
                                    mat = mat,
                                    dir = dir,
                                    rhs = rhs,
                                    bounds = bounds,
                                    max = FALSE,
                                    ...)

  if(denom_lp$optimum + constant_denominator <= 0){
    warning("The minimum value of the denominator is not positive. Returning NaN.")
    return(list(optimum=NaN))
  }

  #We use the Charnes Cooper transformation to convert the linear fractional program into a linear program
  #We add a new optimization variable t>=0
  # The new objective is obj_numerator * (original optimization vector) / obj_denominator + constant_numerator * t
  # The new constraints are :
  # mat * (original optimization_vector) (<=/ == / >=) rhs * t
  # obj_denominator * (original optimization_vector) + constant_denominator * t = 1
  #
  obj_flp <- c(obj_numerator, constant_numerator)

  #Add -rhs as a column to A, then convert rhs to zero
  mat_flp <- base::cbind(mat, -rhs)
  rhs_flp <- rep(0, NROW(rhs))

  #Add the constraint obj_denominator * (orignal optimization vector) + constant_denominator * t = 1
  mat_flp <- rbind(mat_flp,
                   c(obj_denominator, constant_denominator))

  #Update the rhs with this new constraint
  rhs_flp <- c(rhs_flp, 1)

  #Update dir with this new constraint
  dir_flp <- c(dir, "==")

  #Since t>=0 and by default bounds uses bounds of [0,Inf) we don't need to update bounds
  bounds_flp <- bounds

  #Solve the converted fraciton linear program
  fractional_lp_solution <- Rglpk::Rglpk_solve_LP(obj = obj_flp,
                                       mat = mat_flp,
                                       dir = dir_flp,
                                       rhs = rhs_flp,
                                       bounds = bounds_flp,
                                       max = max,
                                       ...)

  #Convert the fractional_lp_solution back to the original by dividing by t (the last element of the optimization vector)
  original_lp_solution <- fractional_lp_solution$solution[-length(fractional_lp_solution$solution)]/fractional_lp_solution$solution[length(fractional_lp_solution$solution)]

  #Return the results
  return(list(optimum = fractional_lp_solution$optimum,
              solution = original_lp_solution,
              status = fractional_lp_solution$status,
              transformed_lp_solution = fractional_lp_solution))
}


compute_max_p_difference <- function(dvec, mdf, yvec, wvec=NULL,...){

  compute_max_p_difference_helper <- function(mvalue){

    #Find all the rows of mdf that equal mvalue
    mindex <- base::apply(mdf, 1, function(x){ base::all(x == mvalue) })

    #Compute density of y among units with m=mvalue and d=1
      #Check if there are any units with m=mvalue and d=1
      #If not, return 0; otherwise, return density
    if(sum(mindex & dvec == 1) == 0){
      dens_y_1 <- function(ygrid){return(0*ygrid)}
    }else{
      dens_y_1 <-
      get_density_fn(x = yvec[mindex & dvec == 1],
                     weights = wvec[mindex & dvec == 1]/sum(wvec[mindex & dvec == 1]),
                     ...)
    }

    #Compute density of y among units with m=mvalue and d=0
      #Check if there are any units with m=mvalue and d=1
      #If not, return 0; otherwise, return density
    if(sum(mindex & dvec == 0) == 0){
      dens_y_0 <- function(ygrid){return(0*ygrid)}
    }else{
      dens_y_0 <-
        get_density_fn(x = yvec[mindex & dvec == 0],
                       weights = wvec[mindex & dvec == 0]/sum(wvec[mindex & dvec == 0]),
                       ...)
    }

    #Compute probability of M=mvalue among D=1 units
    p_m_1 <- stats::weighted.mean( x = mindex[dvec == 1],
                                   weights = wvec)

    p_m_0 <- stats::weighted.mean( x = mindex[dvec == 0],
                                   weights = wvec)

    #Compute integral of max{p_m_1*dens_y_1, p_m_0*dens_y_0} over y
    ygrid <- seq(from = base::min(yvec) - 3* stats::sd(yvec),
                 to = base::max(yvec) + 3* stats::sd(yvec) ,
                 length.out = 10000)

    positive_part <- function(y){ base::pmax(y,0) }

    max_p_diff <- base::sum( positive_part( p_m_1*dens_y_1(ygrid) - p_m_0*dens_y_0(ygrid) ) ) * base::diff(ygrid)[1]

    return(max_p_diff)
  }

  #Compute the max p difference over all unique values of m
  mvalues <- base::unique(mdf)

  max_p_diffs <- base::apply(mvalues, 1, compute_max_p_difference_helper)

  return(list(mvalues = mvalues,
             max_p_diffs = max_p_diffs))
}

compute_partial_densities_and_shares <-
  function(df, d, m, y, w= NULL,...){
    yvec <- df[[y]]
    dvec <- df[[d]]
    mvec <- df[[m]]

    if(is.null(w)){
      wvec <- rep(1, NROW(df))
    }else{
      wvec <- df[[w]]
    }


    frac_compliers <- stats::weighted.mean( mvec[dvec == 1], w = wvec[dvec == 1] ) -
      stats::weighted.mean( mvec[dvec == 0], w = wvec[dvec == 0])
    frac_ats <- stats::weighted.mean( mvec[dvec == 0], w= wvec[dvec == 0] )
    theta_ats <- frac_ats / (frac_compliers + frac_ats) #fraction among Cs/ATs

    ats_untreated_index <- (dvec == 0) & (mvec == 1) #these are ATs when untreated
    ats_treated_index <- (dvec == 1) & (mvec == 1) #these are ATs or Cs

    y_ats_treated <- yvec[ats_treated_index]
    y_ats_untreated <- yvec[ats_untreated_index]

    w_ats_treated <- wvec[ats_treated_index]
    w_ats_untreated <- wvec[ats_untreated_index]

    #The density function doesn't normalize weights, so normalize these
    w_ats_treated <- w_ats_treated/sum(w_ats_treated)
    w_ats_untreated <- w_ats_untreated/sum(w_ats_untreated)

    dens_y_ats_treated <- get_density_fn(x = y_ats_treated, weights = w_ats_treated, ...)
    dens_y_ats_untreated <- get_density_fn(x = y_ats_untreated, weights = w_ats_untreated, ...)

    f_partial11 <- function(y){ (frac_ats + frac_compliers) * dens_y_ats_treated(y) }
    f_partial01 <- function(y){ frac_ats  * dens_y_ats_untreated(y) }

    resultsList <-
      list(frac_compliers = frac_compliers,
           frac_ats = frac_ats,
           theta_ats = theta_ats,
           f_partial11 = f_partial11,
           f_partial01 = f_partial01)

    return(resultsList)

  }


get_density_fn <- function(x,...){

  dens <- stats::density(x=x,...)
  #Create a function that returns the density at any point
  # Return 0 if outside of range of dens$x
  # Otherwise, return the density at the largest pt in dens$x below y
  dens_fn <- function(y){
    index <- base::findInterval(y, dens$x)
    y_plus_zeros <- c(0,dens$y,0)
    return(y_plus_zeros[index+1])
  }

  return(dens_fn)
}

TV_distance_fn <- function(y1,
                           y2,
                           method = "density",
                           numbins = ceiling(length(y1)/20)){

  if(method == "bins"){
    binned_y <- dplyr::ntile( c(y1,y2), numbins)

    binned_y1 <- binned_y[(1:length(y1))]
    binned_y2 <- binned_y[-(1:length(y1))]

    y1_freq_table <- data.frame(y1 = binned_y1) %>%
      dplyr::group_by(y1) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dens1  = n/base::sum(n) )

    y2_freq_table <- data.frame(y2 = binned_y2) %>%
      dplyr::group_by(y2) %>%
      dplyr::count() %>%
      dplyr::ungroup() %>%
      dplyr::mutate(dens2  = n/sum(n) )

    joint_freq_table <- dplyr::full_join(y1_freq_table, y2_freq_table,
                                         by = c("y1" = "y2"))


    joint_freq_table <- joint_freq_table %>%
      dplyr::mutate(dens1 = base::ifelse(is.na(dens1), 0 , dens1),
                    dens2 = base::ifelse(is.na(dens2), 0 , dens2),)

    TV <- 1/2 * base::sum(base::abs(joint_freq_table$dens1 - joint_freq_table$dens2))
    return(TV)
  }
  if(method == "density"){
    dens1 <- stats::density(x = y1)
    dens2 <- stats::density(x = y2)

    dens_fn <- function(dens,y){
      # if(y < min(dens$x) | y> max(dens$x)){
      #   #If Y is out of range, return 0
      #   return(0)
      # }else{
      #   #If Y is in range, find the closest value in the grid returned by density
      #     #This is max value returned by dens such that dens <= y
      #   return(dens[max(which(dens$x <= y))])
      # }

      index <- base::findInterval(y, dens$x)
      y_plus_zeros <- c(0,dens$y,0)
      return(y_plus_zeros[index+1])
    }
    dens1_fn <- function(y){dens_fn(dens1,y)}
    dens2_fn <- function(y){dens_fn(dens2,y)}

    ygrid <- base::seq(from = min(c(dens1$x,dens2$x)),
                       to = max(c(dens1$x,dens2$x)),
                       length.out = 10000)

    #Compute numeric integral of 1/2 |dens1_fn - dens2_fn|
    TV <- 1/2* base::sum( base::abs( dens1_fn(ygrid) - dens2_fn(ygrid) ) * base::diff(ygrid)[1] )
    return(TV)
  }
  else{base::stop("method must be 'density' or 'bins'")}
}



# ##Test TV function
# y1 <- rnorm(10^2, mean = 0, sd = 1)
# y2 <- rnorm(10^2, mean = 1, sd  = 1)
#
# TV_distance_fn(y1,y2)
# TV_distance_fn(y1,y2, method = "bins")
#
# integrate(f = function(.x){ 1/2* abs( dnorm(.x) - dnorm(.x-1) ) } ,
#           lower = -5, upper =5)