#' @title Create plots of partial densities
#' @description Plots f_{Y,M=1 | D=1} and f_{Y,M=1 | D=0}
#' @param df A data frame
#' @param d Name of the treatment variable in the df
#' @param m Name of the mediator variable
#' @param y Name of the outcome variable
#' @param numGridPoints Number of points used in grid for graphing
#' @return A ggplot object showing partial densities
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 xlab
#' @importFrom ggplot2 ylab

#' @export

partial_density_plot <- function(df, d, m, y, numGridPoints = 10000){

  yvec <- df[[y]]

  partial_densities_and_shares <- compute_partial_densities_and_shares(df = df,
                                                                       d = d,
                                                                       m = m,
                                                                       y = y,
                                                                       n = numGridPoints)

  f_partial11 <- partial_densities_and_shares$f_partial11
  f_partial01 <- partial_densities_and_shares$f_partial01

  ygrid <- seq(from = base::min(yvec) - 1* stats::sd(yvec),
               to = base::max(yvec) + 1* stats::sd(yvec) ,
               length.out = numGridPoints)

  partial11_grid <- f_partial11(ygrid)
  partial01_grid <- f_partial01(ygrid)

  partial_density_df <-
  dplyr::bind_rows(
    data.frame(y = ygrid, dens = partial11_grid, Partial.Density = "f(Y,M=1|D=1)"),
    data.frame(y = ygrid, dens = partial01_grid, Partial.Density = "f(Y,M=1|D=0)")
    )

  partial_density_plot <-
  partial_density_df %>%
    ggplot(aes(x = y,
               y = dens,
               color = Partial.Density,
               linetype = Partial.Density)) +
    geom_line() +
    xlab("Y") +
    ylab("Partial Density")

  return(partial_density_plot)
}
