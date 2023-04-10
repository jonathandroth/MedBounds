#' @title Bounds on fraction of always-takers effected
#'@description Computes bounds on the TV distance btwn Y(1,1) and Y(1,0) for always takers
#' This is a lower bound on the fraction of ATs for whom there is a direct effect of D on Y
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable
#' @param method Either "density", to use a kernel density to compute TV, or "bins", to use a discrete approximation
#' @export

compute_tv_ats <- function(df, d, m, y, method = "density"){

  yvec <- df[[y]]
  dvec <- df[[d]]
  mvec <- df[[m]]

  ats_untreated_index <- (dvec == 0) & (mvec == 1)
  ats_treated_index <- (dvec == 1) & (mvec == 1) #these are ATs or Cs

  y_ats_treated <- yvec[ats_treated_index]
  y_ats_untreated <- yvec[ats_untreated_index]

  frac_compliers <- base::mean( mvec[dvec == 1] ) - base::mean( mvec[dvec ==0] )
  frac_ats <- base::mean( mvec[ dvec == 0 ] )
  theta_ats <- frac_ats / (frac_compliers + frac_ats) #fraction among Cs/ATs

  dens_y_ats_treated <- get_density_fn(y_ats_treated)
  dens_y_ats_untreated <- get_density_fn(y_ats_untreated)

  f_partial11 <- function(y){ (frac_ats + frac_compliers) * dens_y_ats_treated(y) }
  f_partial01 <- function(y){ frac_ats  * dens_y_ats_untreated(y) }

  positive_part <- function(y){ base::pmax(y,0) }

  #Lower bound is 1/theta_{AT} \int positive_part( f_partial01 - f_partial11 )
  #Implement this numerically over a grid of 10,000 points

  ygrid <- seq(from = base::min(yvec) - 3* stats::sd(yvec),
               to = base::max(yvec) + 3* stats::sd(yvec) ,
               length.out = 10000)

  TV_lb <-
    1/theta_ats *
  base::sum( positive_part( f_partial01(ygrid) - f_partial11(ygrid) ) ) * base::diff(ygrid)[1]


  # ##Old approach using non-sharp bounds
  # #Compute TV distance between outcome when D=1,M=1 versus D=0,M=1
  # TV_11_vs_10 <- TV_distance_fn( y_ats_treated,
  #                                y_ats_untreated,
  #                                method = method)
  #
  # #Lower bound on TV of Y(1,1) vs Y(0,1) for ATs is:
  # # 1/theta_ats * (TV_11_vs_10 - (1-theta_ats))
  # TV_lb <- max(1/theta_ats * (TV_11_vs_10 - (1-theta_ats)), 0)
  #
  return(TV_lb)
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