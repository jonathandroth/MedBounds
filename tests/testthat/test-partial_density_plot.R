test_that("Partial density plot returns a ggplot", {
  expect_equal(class(partial_density_plot(kerwin_data,
                                          d= "treated",
                                          m = "primarily_leblango",
                                          y = "EL_EGRA_PCA_Index"))[2], "ggplot")
})
