#' @title Hypothesis test for the sharp null
#' @description This function tests the sharp null of Y(1,m) = Y(0,m). The
#'   outcome and mediator are both assumed to take finitely many different
#'   values. The inference is via applying Cox and Shi (2023)
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable, which is assumed to take a discrete
#'   support
#' @param ordering A list with length equal to the cardinality of the support of
#'   the mediator variable. The name of each element corresponds to a point in
#'   the support, and each element is a vector that collects all m values that
#'   are less than or equal to this point. If ordering = NULL, the standard
#'   ordering is used.
#' @param B Bootstrap size, default is zero
#' @param cluster Cluster for bootstrap
#' @param weight.matrix Weight matrix used to implement FSST. Possible options
#'   are "diag", "avar", "identity." Defaults is "diag" as in FSST.
#'  @param num_Ybins (Optional) If specified, Y is discretized into the given number of bins (if num_Ybins is larger than the number of unique values of Y, no changes are made)
#'  @param refinement (Optional) If TRUE, use the refined Cox & Shi test (rCC rather than CC). Default is FALSE.
#' @export
test_sharp_null_coxandshi_binary_m <- function(df,
                                         d,
                                         m,
                                         y,
                                         ordering = NULL,
                                         B = 500,
                                         cluster = NULL,
                                         weight.matrix = "diag",
                                         ats_only = F,
                                         alpha = 0.05,
                                         kappa = alpha/10,
                                         use_hybrid = T,
                                         num_Ybins = NULL,
                                         analytic_variance = FALSE,
                                         refinement = FALSE,
                                         print_both_var = FALSE){


  df <- remove_missing_from_df(df = df,
                               d = d,
                               m = m,
                               y = y)

  yvec <- df[[y]]

  if(!is.null(num_Ybins)){
    yvec <- discretize_y(yvec = yvec, numBins = num_Ybins)
    df[[y]] <- yvec
  }

  if (is.null(cluster)) {
    clustervec <- 1:length(yvec)
  } else {
    clustervec <- df[[cluster]]
  }

  dvec <- df[[d]]
  mvec <- df[[m]]


  yvalues <- sort(unique(yvec))
  mvalues <- unique(mvec)
  my_values <- purrr::cross_df(list(m=mvalues,y=yvalues)) %>%
    dplyr::arrange(m,y) %>%
    dplyr::select(y,m)

  get_beta.obs <- function(yvec, dvec, mvec) {
    #Get partial density for Y,M=1|D=1
    p_y1_1 <- purrr::map_dbl(.x = 1:length(yvalues),
                             .f = ~mean(yvec[dvec == 1] == yvalues[.x]
                                        & mvec[dvec == 1] == 1 ))

    #Get partial density for Y,M=1|D=0
    p_y1_0 <- purrr::map_dbl(.x = 1:length(yvalues),
                             .f = ~mean(yvec[dvec == 0] == yvalues[.x]
                                        & mvec[dvec == 0] == 1 ))

    if(!ats_only){
      #Get partial density for Y,M=0|D=1
      p_y0_1 <- purrr::map_dbl(.x = 1:length(yvalues),
                               .f = ~mean(yvec[dvec == 1] == yvalues[.x]
                                          & mvec[dvec == 1] == 0 ))

      #Get partial density for Y,M=0|D=0
      p_y0_0 <- purrr::map_dbl(.x = 1:length(yvalues),
                               .f = ~mean(yvec[dvec == 0] == yvalues[.x]
                                          & mvec[dvec == 0] == 0 ))

    }
    #We return differences in partial densities that should be positive
    if(ats_only){
      beta.obs <- c(p_y1_1 - p_y1_0)
    }else{
      beta.obs <- c(p_y0_0 - p_y0_1,
                    p_y1_1 - p_y1_0)

    }
    return(beta.obs)
  }

  #Bootstrap the betas
  beta.obs_list <- compute_bootstrap_draws_clustered(f =
                                                       function(df,d,y,m,...){get_beta.obs(
                                                         df[[y]],
                                                         df[[d]],
                                                         df[[m]])},
                                                     df = df,
                                                     d = d,
                                                     m = m,
                                                     y = y,
                                                     cluster = cluster,
                                                     numdraws = B,
                                                     return_df = F)

  # Get variance matrix of the beta.obs bootsraps
  if (analytic_variance) {
    #Calculate the analytic variance
    sigma.obs <- analytic_variance(yvec = yvec,
                                   dvec = dvec,
                                   mvec = mvec,
                                   my_values = my_values,
                                   inequalities_only = TRUE,
                                   clustervec = clustervec,
                                   exploit_binary_m = TRUE)


    if (print_both_var) {
      sigma.obs_boot <- stats::cov(base::Reduce(base::rbind,
                                                beta.obs_list))
      print(sigma.obs)
      print(sigma.obs_boot)
    }
  } else {
    sigma.obs <- stats::cov(base::Reduce(base::rbind,
                                       beta.obs_list))
  }


  #Get the beta.obs using actual data
  beta.obs <- get_beta.obs(yvec, dvec, mvec)

  #Run Cox and Shi test
  coxandshi <- cox_shi_nonuisance(Y = -beta.obs,
                                  sigma = sigma.obs,
                                  alpha = alpha,
                                  refinement = refinement)
  return(coxandshi)
}
