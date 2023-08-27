#' @description This function processes the inputted dataframe to remove missing values
#' @param df A data.frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable
#' @param w (Optional) Name of weighting variable. If null, equal weights are used

remove_missing_from_df <- function(df, d, m, y, w = NULL){

  missing_d <- is.na(df[[d]])
  missing_m <- pmax(is.na(df[[m]])) #allow for vector-valued m
  missing_y <- is.na(df[[y]])

  df <- df[!(missing_d | missing_m | missing_y), ]

  if(!is.null(w)){
    missing_w <- is.na(df[[w]])
    df <- df[!missing_w, ]
  }

  return(df)
}
