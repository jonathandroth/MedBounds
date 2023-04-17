#' @title Bounds on direct for always-takers
#' @description Computes bounds on E[Y(1,1) - Y(1,0) | G = AT], the average direct effect for ATs for whom there is a direct effect of D on Y
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable
#' @importFrom "stats" "quantile"
#' @export

compute_bounds_ats <- function(df, d, m, y, w = NULL){

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

  ### ATs Y(1,1) outcome is partially-identified ###
  ats_treated_index <- (dvec == 1) & (mvec == 1) #these are ATs or Cs

  #Compute fraction of ATs/Cs
  frac_compliers <- stats::weighted.mean( mvec[dvec == 1], w =  wvec[dvec == 1] ) -
    stats::weighted.mean( mvec[dvec == 0], w = wvec[dvec == 0] )
  frac_ats <- stats::weighted.mean( mvec[ dvec == 0 ], w = wvec[ dvec == 0] )
  theta_ats <- frac_ats / (frac_compliers + frac_ats) #fraction among Cs/ATs

  # Get ub on Y(1,1) by trimming to the upper theta_at's fraction
  quantile_ub <- Hmisc::wtd.quantile(x = yvec[ats_treated_index],
                                     probs = 1-theta_ats,
                                     weights = wvec[ats_treated_index],
                                     normwt = T)

  ub_treated_mean <- stats::weighted.mean( x = yvec[ ats_treated_index & yvec >= quantile_ub  ],
                                           w = wvec[ ats_treated_index & yvec >= quantile_ub  ])

  # Get lb on Y(1,1) by trimming to the lower theta_at's fraction
  quantile_lb <- Hmisc::wtd.quantile(x = yvec[ats_treated_index],
                                     probs = theta_ats,
                                     weights = wvec[ats_treated_index],
                                     normwt = T)
  lb_treated_mean <- stats::weighted.mean( x = yvec[ ats_treated_index & yvec <= quantile_lb  ],
                                           w = wvec[ ats_treated_index & yvec <= quantile_lb  ])


  #Get bounds on treatment effect by taking up/lbs for treated minus control mean
  ub <- ub_treated_mean - ats_untreated_mean
  lb <- lb_treated_mean - ats_untreated_mean
  return( data.frame(lb = lb, ub = ub) )
}
