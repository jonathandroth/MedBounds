#' @title Bounds on fraction of always-takers affected with a multi-valued discrete or multi-dimensional discrete m
#'@description Computes  a lower bound on the TV distance btwn Y(1,m) and Y(0,m) for always takers
#' who have M(1)=M(0)=m if at_group=m. This is a lower bound on the fraction of these ATs for whom there is a direct effect of D on Y
#' If at_group=NULL, it computes a lower bound on the weighted average across all ATs, weighted by their proportion in the population
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Vector of the mediator variable
#' @param y Name of the outcome variable
#' @param w (Optional) Name of weighting variable. If null, equal weights are used
#' @param at_group (Optional) Value of m specifying which always-takers to compute lower bounds of TV for.
#' If at_group is specified, then we compute a lower bound on TV between Y(1,at_group) and Y(0,at_group) for
#' ATs who have M(1)=M(0)=at_group. If at_group is null (the default), we compute a lower bound on
#' the weighted average of TV across all always-takers, with weights proportional to shares in population
#' @param continuous_Y (Optional) Whether Y should be treated as continuous, in which case kernel density is used, or discrete. Default is TRUE.
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
#' @param max_defier_share (Optional) If specified, up to max_defier_share fraction of the population can be defiers (of any type). Default is zero
#' @param allow_min_defiers (Optional) If the bound on defiers (max_defier_share) is inconsistent with the data, proceed by allowing the minimum number of defiers compatible with the data. Otherwise, throw an error. Default is FALSE
#' @param return_min_defiers (Optional) If true, the function returns the minimum number of defiers consistent with the data rather than the TV bounds
#' @export

compute_tv_ats_multiple_m <- function(df,
                                      d,
                                      m,
                                      y,
                                      at_group = NULL,
                                      w = NULL,
                                      continuous_Y = base::ifelse(is.null(num_Ybins),TRUE,FALSE),
                                      num_Ybins = NULL,
                                      max_defier_share = 0,
                                      allow_min_defiers = FALSE,
                                      return_min_defiers = FALSE){

  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y,
                               w = w)


  yvec <- df[[y]]

  if(!is.null(num_Ybins)){
    yvec <- discretize_y(yvec = yvec, numBins = num_Ybins)
    df[[y]] <- yvec
  }

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
                                               wvec=wvec,
                                               continuous_Y = continuous_Y)
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

    defier_types <- 1-monotonic_types

    #Restrict to the types satisfying monotonicity
    #m0_types <- m0_types[monotonic_types,]
    #m1_types <- m1_types[monotonic_types,]


    #Compute the marginal distribution of M among D=1
    p_m_1_fn <- function(mvalue){
      mdf1 <- mdf[dvec ==1,]
      M_equals_mvalue <- sapply(1:NROW(mdf1), function(i){all(mdf1[i,]==mvalue)} )
      return(stats::weighted.mean( x = M_equals_mvalue,
                                   w = wvec[dvec == 1]/sum(wvec[dvec == 1])))
    }

    p_m_1 <- base::apply(mvalues,1, p_m_1_fn)

    #Compute the marginal distribution of M among D=0
    p_m_0_fn <- function(mvalue){
      mdf0 <- mdf[dvec ==0,]
      M_equals_mvalue <- sapply(1:NROW(mdf0), function(i){all(mdf0[i,]==mvalue)} )
      return(stats::weighted.mean( x = M_equals_mvalue,
                                   w = wvec[dvec == 0]/sum(wvec[dvec == 0])))
    }

    p_m_0 <- base::apply(mvalues,1, p_m_0_fn)

    ## We now calculate lower bounds on TV using either the group provided, or averaging across groups proportionally to share
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

    #Create a matrix (really, vector) that bounds the total defiers share
    defiers_constraints_matrix <- c(defier_types, rep(0,NROW(mvalues)))


    ## We now check feasibility of the program by computing the minimal defier share
    # consistent with the constraints
    feasibility_lp <-
      Rglpk::Rglpk_solve_LP(obj = defiers_constraints_matrix, #obj is sum of defiers shares
                            mat = rbind(m1_marginals_constraints_matrix,
                                        m0_marginals_constraints_matrix),
                            rhs = c(p_m_1,p_m_0),
                            dir = rep("==", NROW(m1_marginals_constraints_matrix)*2),
                            max = FALSE
      )

    if(feasibility_lp$status == 1){
      warning("Error in checking feasibility. Proceed with caution")
    }else{

      min_defier_share <- feasibility_lp$optimum

      #Return the min defier share if specified
      if(return_min_defiers){return(min_defier_share)}

      if(min_defier_share > max_defier_share){
        if(allow_min_defiers){
          max_defier_share <- min_defier_share + 10^(-6) #update max defier share to the optimum plus a small tolerance
          warning(paste0("The data is incompatible with the specified max_defier_share.
                           Setting this to the min value compatible with the data:",
                         max_defier_share))
        }else{
          stop(paste0("The data is incompatible with the specified max_defier_share. \n
                        The specified value is ",
                      max_defier_share, " but the min compatible with the data is ", min_defier_share))
        }
      }
    }



    ##Create a constraint corresponding to the constraint that max_p_diffs must be greater than or equal to the sum of all complier types that end up at a given mvalue
    #Create a matrix where each row i is a NROW(m1_types) length vector where the jth element indicates if the jth row of m1_types equals the ith row of mvalues and the jth row of m1_types is not exactly equal to the jth row of m0_types
    complier_constraints_matrix <-
      Reduce(rbind,
             purrr::map(.x = 1:NROW(mvalues),
                        .f = function(m_index){ purrr::map_lgl(.x = 1:NROW(m1_types),
                                                               .f = ~MedBounds:::row_equals(m1_types[.x,], mvalues[m_index,] ) &
                                                                     !MedBounds:::row_equals(m0_types[.x,], mvalues[m_index,] ))  })
            )
        # base::lapply(X = 1:NROW(mvalues),
        #              FUN = function(i){base::sapply(1:NROW(m1_types), function(s){ base::all(m1_types[s,] == mvalues[i,]) & !base::all(m0_types[s,] == m0_types[i,]) })}))

    #Add an identity matrix to complier_constraints_matrix
      # Thus the end vector corresponds to violation of kth moment, i.e. theta_kk TV_k

    complier_constraints_matrix <- base::cbind(complier_constraints_matrix,
                                               diag(NROW(mvalues)))


    #Combine the constraint matrices
    constraints_matrix <- base::rbind(m1_marginals_constraints_matrix,
                                      m0_marginals_constraints_matrix,
                                      complier_constraints_matrix,
                                      defiers_constraints_matrix)


    #Combine the constants associated with the matrices
    d <- c(p_m_1, p_m_0, max_p_diffs, max_defier_share)

    #Specify the direction of the equalities/inequalities
    dir <- c(rep("==", 2*NROW(m1_marginals_constraints_matrix)),
             rep(">=", NROW(max_p_diffs)),
             "<=")



    #We do not need to provide a bound argument, since by default Rglpk requires variables to be positive
    # The marginal probabilities constraints imply that the the thetas can't be >1
    #(We previously specified bounds here, but they messed up the transformation of the fractional-lp)
    #If you want to provide bounds, these should be added to the matrix of constraints



    ##Specify the objective

    ##Create a function that returns the objective vectors for a single at_group
    ##Our objective is theta_kk TV_k / theta_kk
    compute_obj_vectors <- function(at_group){


      #Find the index i such that the ith row of m0_types and m1_types are equal to at_group
      at_group_index <- which(purrr::map_lgl(.x = 1:NROW(m0_types),
                                             .f = ~MedBounds:::row_equals(m0_types[.x,], at_group) &
                                               MedBounds:::row_equals(m1_types[.x,], at_group)))

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
      #Hence, we compute the objective vector for each at_group and take the average
      obj_numerator <- Reduce("+", base::lapply(X = 1:NROW(mvalues),
                                                FUN = function(i){compute_obj_vectors(at_group = mvalues[i,])$obj_numerator}))
      obj_denominator <- Reduce("+", base::lapply(X = 1:NROW(mvalues),
                                                 FUN = function(i){compute_obj_vectors(at_group = mvalues[i,])$obj_denominator}))
    }
    #Use Rglpk to solve the fraction-linear program to minimize the weighted average of TV
      #We capture warnings if the denom can be zero
    quiet_fractional_LP <- purrr::quietly(MedBounds:::Rglpk_solve_fractional_LP)
    max_violation <-
      quiet_fractional_LP(
        obj_numerator = obj_numerator,
        obj_denominator = obj_denominator,
        mat = constraints_matrix,
        dir = dir,
        rhs = d,
        bounds = NULL,
        max = FALSE)

    if(!is.na(max_violation$warnings[1]) & max_violation$warnings[1] == "The minimum value of the denominator is not positive. Returning NaN."){
      warning("The lower bound for the fraction of always-takers for the chosen at_group is zero. Returning NaN")
      return(NaN)
    }
    #Return the maximal violation
    return(max_violation$result$optimum)
}

#' @title Function that computes whether a row of a df x equals a row of a df y in all elements
row_equals <- function(x,y){ base::all(x == y)}

#' @title Wrapper for Rglpk_solve_lp that implements linear fractional programming
#' @param obj_numerator The objective vector for the numerator of the linear fractional program
#' @param obj_denominator The objective vector for the denominator of the linear fractional program. Note that the denominator must be strictly positive over all feasible values (if not, we return a warning and NaN for the value)
#' @param mat The constraint matrix
#' @param dir The direction of the constraints
#' @param rhs The right hand side of the constraints
#' @param bounds The bounds on the optimization vector. Note: the bounds should only be 0 or Inf, otherwise this will mess up the transformation of the fractional-LP to an LP. If you have other bounds, add these to the constraints in mat
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

  #We first check feasibility
  if(denom_lp$status == 1){
    stop("The LP is not feasible, likely because monotonicity is not satisfied.")
  }

  #Subject to feasibility, we check whether the solution for the denominator is positive
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


compute_max_p_difference <- function(dvec, mdf, yvec, wvec=NULL,
                                     continuous_Y = TRUE,...){

  compute_max_p_difference_helper <- function(mvalue){
      #Find all the rows of mdf that equal mvalue
    mindex <- base::apply(mdf, 1, function(x){ base::all(x == mvalue) })

    #We now compute the partial densities / PMFs, depending on whether Y is cts

    if(continuous_Y == TRUE){
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

    }else{
      yvalues <- unique(yvec[mindex])
      pmf_y_1 <- purrr::map_dbl(.x = 1:length(yvalues),
                                .f = ~stats::weighted.mean(x = yvec[mindex & dvec==1]== yvalues[.x],
                                                          w = wvec[mindex & dvec==1]) )

      pmf_y_0 <- purrr::map_dbl(.x = 1:length(yvalues),
                                .f = ~stats::weighted.mean(x = yvec[mindex & dvec==0]== yvalues[.x],
                                                          w = wvec[mindex & dvec==0]) )


    }

    #Compute probability of M=mvalue among D=1 units
    p_m_1 <- stats::weighted.mean( x = mindex[dvec == 1],
                                   w = wvec[dvec == 1])

    p_m_0 <- stats::weighted.mean( x = mindex[dvec == 0],
                                   w = wvec[dvec == 0])

    if(continuous_Y == TRUE){
    #If continuous Compute integral of max{p_m_1*dens_y_1, p_m_0*dens_y_0} over y
    ygrid <- seq(from = base::min(yvec) - 3* stats::sd(yvec),
                 to = base::max(yvec) + 3* stats::sd(yvec) ,
                 length.out = 10000)

    positive_part <- function(y){ base::pmax(y,0) }

    max_p_diff <- base::sum( positive_part( p_m_1*dens_y_1(ygrid) - p_m_0*dens_y_0(ygrid) ) ) * base::diff(ygrid)[1]
    }else{
      #If discrete, compute positive part of partial PMF diffs
      partial_pmf_1 <- p_m_1 * pmf_y_1
      partial_pmf_0 <- p_m_0 * pmf_y_0

      max_p_diff <- base::sum(base::pmax(partial_pmf_1-partial_pmf_0,0))
    }

    return(max_p_diff)
  }

  #Compute the max p difference over all unique values of m
  mvalues <- base::unique(mdf)

  max_p_diffs <- base::apply(mvalues, 1, compute_max_p_difference_helper)

  return(list(mvalues = mvalues,
             max_p_diffs = max_p_diffs))
}


#' @export
#' @title Finds the minimum number of defiers compatible with the sharp null
#'@description This function finds the minimum value of max_defier_share such that compute_tv_ats_multiple_m returns zero
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Vector of the mediator variable
#' @param y Name of the outcome variable
#' @param w (Optional) Name of weighting variable. If null, equal weights are used
#' @param at_group (Optional) Value of m specifying which always-takers to compute lower bounds of TV for.
#' If at_group is specified, then we compute a lower bound on TV between Y(1,at_group) and Y(0,at_group) for
#' ATs who have M(1)=M(0)=at_group. If at_group is null (the default), we compute a lower bound on
#' the weighted average of TV across all always-takers, with weights proportional to shares in population
#' @param continuous_Y (Optional) Whether Y should be treated as continuous, in which case kernel density is used, or discrete. Default is TRUE.
#' @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)

breakdown_defier_share <- function(df,
                                   d,
                                   m,
                                   y,
                                   at_group = NULL,
                                   w = NULL,
                                   continuous_Y = base::ifelse(is.null(num_Ybins),TRUE,FALSE),
                                   num_Ybins = NULL){

  tol <- 10^(-4)

  #Function for computing the lb as a function of max_defier_share
  lb_fn <- function(max_defier_share){
    min_tv <-
    base::suppressWarnings(
    compute_tv_ats_multiple_m(df = df,
                              d = d,
                              m = m,
                              y = y ,
                              at_group = at_group,
                              w = w,
                              continuous_Y = continuous_Y,
                              num_Ybins = num_Ybins,
                              max_defier_share = max_defier_share))

    #If the lb on ATs is zero, compute_tv_ats_multiple_m returns zero
    #For our purposes, we treat this as a zero
    if(is.nan(min_tv)){min_tv <- 0}
    return(min_tv)
  }

  # Compute the min defier share compatible with the data
  min_compatible_defiers <-     compute_tv_ats_multiple_m(df = df,
                                                          d = d,
                                                          m = m,
                                                          y = y ,
                                                          at_group = at_group,
                                                          w = w,
                                                          continuous_Y = continuous_Y,
                                                          num_Ybins = num_Ybins,
                                                          max_defier_share = max_defier_share,
                                                          return_min_defiers = TRUE)



  #Check if lb is positive under min number of defiers
    #If not, return min number of defiers
  if(lb_fn(min_compatible_defiers) < tol){return(min_compatible_defiers)}

  #Check if lb is positive with max defier share of 1
    #If so, return Inf
  if(lb_fn(1) > tol){return(Inf)}

  #Otherwise, we have a zero somewhere between the min compatible defiers and 1
  # We find it using binary search
  ub <- 1
  lb <- min_compatible_defiers
  while(ub - lb > tol){
    mid <- (ub + lb )/2
    val <- lb_fn(mid)

    if(val > tol){lb <- mid}
    else{ub <- mid}
  }
  return(ub)
}
