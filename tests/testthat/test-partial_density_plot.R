test_that("Partial density plot returns a ggplot", {
  expect_equal(class(partial_density_plot(kerwin_data,
                                          d= "treated",
                                          m = "primarily_leblango",
                                          y = "EL_EGRA_PCA_Index"))[2], "ggplot")
})

test_that("Partial density plot returns a ggplot with plot_nts option", {
  expect_equal(class(partial_density_plot(kerwin_data %>%
                                            dplyr::mutate(primarily_leblango = 1-primarily_leblango,
                                                   treated = 1-treated) , #flip m and d so that we have some NTs
                                          d= "treated",
                                          m = "primarily_leblango",
                                          y = "EL_EGRA_PCA_Index",
                                          plot_nts = T))[2], "ggplot")
})


test_that("Partial density plot returns a ggplot with discrete Y", {
  expect_equal(class(partial_density_plot(kerwin_data %>%
                                            dplyr::mutate(EL_EGRA_PCA_Index = dplyr::ntile(EL_EGRA_PCA_Index,3)), #create a continuous Y ,
                                          d= "treated",
                                          m = "primarily_leblango",
                                          y = "EL_EGRA_PCA_Index",
                                          continuous_Y = FALSE))[2], "ggplot")
})
