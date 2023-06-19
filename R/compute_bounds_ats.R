#' @title Bounds on direct effect for always-takers
#' @description Computes bounds on E[Y(1,1) - Y(1,0) | G = AT], the average direct effect for ATs for whom there is a direct effect of D on Y
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable
#' @param w (Optional) Name of weighting variable. If null, equal weights are used
#' @param c_at_ratio (Optional) specify the ratio of E[Y(1,1) | C]/E[Y(1,1) | AT]. If this is specified, direct effect for ATs is point-identified
#' @param adjust_for_point_mass (Optional) specify whether to use a correction to the bounds that allows for point-mass in the distribution of Y (default TRUE)
#' @importFrom "stats" "quantile"
#' @export

compute_bounds_ats <- function(df,
                               d,
                               m,
                               y,
                               w = NULL,
                               c_at_ratio = NULL,
                               adjust_for_point_mass = TRUE){

  yvec <- df[[y]]
  dvec <- df[[d]]
  mvec <- df[[m]]

  if(is.null(w)){
    wvec <- rep(1, NROW(df))
  }else{
    wvec <- df[[w]]
  }


  #ATs Y(0,1) outcome is point-identified based on observations with D=0,M=1
  ats_untreated_index <- (dvec == 0) & (mvec == 1)
  ats_untreated_mean <- stats::weighted.mean( x = yvec[ats_untreated_index],
                                              w = wvec[ats_untreated_index])



  ##Compute fraction of ATs/Cs
  ats_treated_index <- (dvec == 1) & (mvec == 1) #these are ATs or Cs
  frac_compliers <- stats::weighted.mean( x = mvec[dvec == 1], w = wvec[dvec == 1] ) -
    stats::weighted.mean( x = mvec[dvec == 0], w = wvec[dvec == 0] )
  frac_ats <- stats::weighted.mean( x = mvec[ dvec == 0 ],  w = wvec[ dvec == 0 ]  )
  theta_ats <- frac_ats / (frac_compliers + frac_ats) #fraction among Cs/ATs

  ### ATs Y(1,1) outcome is partially-identified if ratio bound is not specified###

  if(is.null(c_at_ratio)){

  # Get ub on Y(1,1) by trimming to the upper theta_at's fraction
  quantile_ub <- reldist::wtd.quantile(x = yvec[ats_treated_index],
                                       q = 1-theta_ats,
                                       weight = wvec[ats_treated_index])
    # quantile_ub <- Hmisc::wtd.quantile(x = yvec[ats_treated_index],
    #                                    probs = 1-theta_ats,
    #                                    weights = wvec[ats_treated_index],
    #                                    normwt = TRUE)

  if(adjust_for_point_mass){
    #Quantile function uses interpolation. We want the inf version if we're going to correct
    quantile_ub <- min(yvec[ats_treated_index &  yvec >= quantile_ub])
  }

  ub_treated_mean <- stats::weighted.mean( yvec[ ats_treated_index & yvec >= quantile_ub  ],
                                           wvec[ ats_treated_index & yvec >= quantile_ub  ])


  #Adjust the trimmed mean to allow for point-mass at quantile_ub
  if(adjust_for_point_mass == TRUE){
    p <- stats::weighted.mean( x = (yvec[ats_treated_index] < quantile_ub ),
                               w = wvec[ats_treated_index] )

    thetatilde <- 1 - (1-theta_ats - p)/(1-p)
    ub_treated_mean <- 1/thetatilde * ub_treated_mean - (1-thetatilde) / thetatilde * quantile_ub
  }

  # Get lb on Y(1,1) by trimming to the lower theta_at's fraction
  quantile_lb <- reldist::wtd.quantile(x = yvec[ats_treated_index],
                                       q = theta_ats,
                                       weight = wvec[ats_treated_index])

  # quantile_lb <- Hmisc::wtd.quantile(x = yvec[ats_treated_index],
  #                                    probs = theta_ats,
  #                                    weights = wvec[ats_treated_index],
  #                                    normwt = TRUE)

  if(adjust_for_point_mass){
  #Quantile function uses interpolation. We want the inf version if we're going to correct
  quantile_lb <- min(yvec[ats_treated_index &  yvec >= quantile_lb])
  }

  lb_treated_mean <- stats::weighted.mean( x = yvec[ ats_treated_index & yvec <= quantile_lb  ],
                                           w = wvec[ ats_treated_index & yvec <= quantile_lb  ])


  #Adjust the trimmed mean to allow for point-mass at quantile_ub
  if(adjust_for_point_mass == TRUE){
    p <- stats::weighted.mean( x = (yvec[ats_treated_index] <= quantile_lb ),
                               w = wvec[ats_treated_index] )

    thetatilde <- theta_ats / p
    lb_treated_mean <- 1/thetatilde * lb_treated_mean - (1-thetatilde) / thetatilde * quantile_lb
  }

  #Get bounds on treatment effect by taking up/lbs for treated minus control mean
  ub <- ub_treated_mean - ats_untreated_mean
  lb <- lb_treated_mean - ats_untreated_mean
  return( data.frame(lb = lb, ub = ub) )

  }else{
    mean_for_d1m1 <- stats::weighted.mean( x = yvec[ats_treated_index],
                                           w = wvec[ats_treated_index])
    ats_treated_mean <- 1/(theta_ats + c_at_ratio * (1- theta_ats)) * mean_for_d1m1
    ate_ats <- ats_treated_mean - ats_untreated_mean
    return(data.frame(ate_ats = ate_ats))
  }
}

compute_posterior_draws <- function(f, df, d, m, y, w = NULL, numdraws = 100){
  n <- NROW(df)

  compute_posterior_helper <- function(seed){
    set.seed(seed)
    gammaDraws <- stats::rgamma(n = n, shape = 1) #draw indep gammas
    dirichletDraws <- gammaDraws / sum(gammaDraws) #normalize to get Dirichlet
    df$dirichletDraws <- dirichletDraws
    return(f(df, d, m, y, w = "dirichletDraws"))
  }

  posteriorDraws <- purrr::map_dfr(.x = 1:numdraws, .f = compute_posterior_helper)
  return(posteriorDraws)
}
