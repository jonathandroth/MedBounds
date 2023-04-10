test_that("Lee bounds and ours coincide", {
  library(leebounds)
  leeresults <-
    leebounds::leebounds(MedBounds::kerwin_data %>%
                dplyr::rename(outcome = EL_EGRA_PCA_Index,
                       selection = primarily_leblango,
                       treat = treated))

  bounds_ats <-
  compute_bounds_ats(df = kerwin_data,
                     d= "treated",
                     m = "primarily_leblango",
                     y = "EL_EGRA_PCA_Index")
  expect_equal(bounds_ats$lb, leeresults$lower_bound,
               tolerance = 0.001)

  expect_equal(bounds_ats$ub, leeresults$upper_bound,
               tolerance = 0.001)
})
