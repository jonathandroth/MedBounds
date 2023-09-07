#' @title Function for discretizing yvec.
#' @param yvec A vector containing y values (assumed non-missing)
#' @param numBins The target number of bins into which to discretize y
#' @returns A discretized version of yvec. If yvec has fewer than numBins unique
#'  values, then yvec is returned. Otherwise, create cutpoints for y using the
#'  1/numBins, 2/numBins.... quantiles of y. If two cutpoints are the same (owing to point mass)
#'  we group those together into one bin, so the number of bins returns may be less than numBins

discretize_y <- function(yvec, numBins){
  yvalues <- base::unique(yvec)

  if(length(yvalues) <= numBins){
    return(yvec)
  }

  cutpoints <- stats::quantile(x = yvec,
                               prob =
                                 seq(from=1, to = numBins-1, by = 1)/numBins)

  #Add -Inf and +Inf to the cutpoints so that all values are include
  # And restrict to unique values of cutpoints
  cutpoints <- c(-Inf, unique(cutpoints), Inf)

  yvec_discretized <- base::cut(x = yvec,
                                breaks = cutpoints)

  return(yvec_discretized)
}
