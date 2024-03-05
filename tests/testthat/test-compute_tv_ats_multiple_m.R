#Write unit tests for the compute_tv_ats_multiple_m function

#compute_tv_ats_multiple_m should match the output of compute_tv_ats when a single binary M is provided

test_that("compute_tv_ats_multiple_m matches compute_tv_ats when a binary M is provided", {
  #Use the fn for the binary M case
  tv_binary <-
  MedBounds::compute_tv_ats(df = MedBounds::kerwin_data,
                              d= "treated",
                              m = "primarily_leblango",
                              y = "EL_EGRA_PCA_Index")

  tv_multiple <-
    MedBounds::compute_tv_ats_multiple_m(df = MedBounds::kerwin_data,
                              d= "treated",
                              m = "primarily_leblango",
                              y = "EL_EGRA_PCA_Index",
                              at_group = c(TRUE,TRUE))

  #Check that the two are equal
  expect_equal(tv_binary, tv_multiple, tolerance = 0.01)
})


#compute_tv_ats_multiple_m should agree with the analytic formulas when M is an ordered scalar
test_that("compute_tv_ats_multiple_m should agree with the analytic formulas when M is an ordered scalar",{
  #Create an M with three values
  df <- MedBounds::kerwin_data
  df$m <- dplyr::ntile(df$primarily_leblango,3)

  #Calculate TV_{22} using analytic formulas
  p_m2_1 <- mean(df$m[df$treated == 1] == 2 )
  p_mgte2_1 <- mean(df$m[df$treated == 1] >= 2 )
  p_mgte2_0 <- mean(df$m[df$treated == 0] >= 2 )

  theta_kk <- p_m2_1 - (p_mgte2_1 - p_mgte2_0)

  max_p_diffs_list <- compute_max_p_difference(dvec = df$treated,
                                               mdf = df[,"m"],
                                               yvec = df$EL_EGRA_PCA_Index,
                                               wvec = rep(1/NROW(df),NROW(df)))

  max_p_diff <- max_p_diffs_list$max_p_diffs[max_p_diffs_list$mvalues == 2]

  analytical_tv2 <- (max_p_diff - (p_mgte2_1-p_mgte2_0))/ theta_kk

  tv_multiple <-
    MedBounds::compute_tv_ats_multiple_m(df = df,
                                         d= "treated",
                                         m = "m",
                                         y = "EL_EGRA_PCA_Index",
                                         at_group = 2)

  expect_equal(analytical_tv2, tv_multiple, tolerance = 0.01)

})


#Test the fractional linear programming wrapper that we wrote
test_that("Fractional linear programming works on simple examples", {

  frac_lp_result <-
  MedBounds:::Rglpk_solve_fractional_LP(obj_numerator = c(1,0),
                            obj_denominator = c(0,1),
                            constant_numerator = 0,
                            constant_denominator = 0,
                            mat = matrix(c(1,0,0,-1),byrow=T,nrow=2),
                            rhs = c(3,-1/2),
                            dir=c("<=","<="),
                            max = TRUE, bounds = NULL
  )

  expect_equal(frac_lp_result$optimum,6, tolerance = 0.01)
  expect_equal(max(abs(frac_lp_result$solution - c(3,1/2))),0, tolerance = 0.01)

  frac_lp_result <-
    MedBounds:::Rglpk_solve_fractional_LP(obj_numerator = c(1,0),
                                          obj_denominator = c(0,1),
                                          constant_numerator = 0,
                                          constant_denominator = 0,
                                          mat = diag(2),
                                          rhs = c(3,1),
                                          dir=c("<=",">="),
                                          max = TRUE, bounds = NULL
    )

  expect_equal(frac_lp_result$optimum, 3, tolerance = 0.01)

  frac_lp_result <- MedBounds:::Rglpk_solve_fractional_LP(obj_numerator = c(-2,1),
                                                          obj_denominator = c(1,3),
                                                          constant_numerator = 2,
                                                          constant_denominator = 4,
                                                          mat = matrix(c(-1,1,2,1),byrow=T,nrow=2),
                                                          rhs = c(4,14),
                                                          dir=c("<=","<="),
                                                          bounds = list(lower=list(ind = c(1,2), val = c(0,0)),
                                                                        upper=list(ind = c(1,2), val = c(Inf,6))),
                                                          max = F)
  expect_equal(frac_lp_result$optimum, -12/11, tolerance = 0.01)
  expect_equal( max(abs(frac_lp_result$solution- c(7,0))), 0, tolerance = 0.01)
})

test_that("get error when monotonicity is violated",{
  expect_error(compute_tv_ats_multiple_m(
    kerwin_data %>% dplyr::mutate(minus_treated = 1-treated),
    d = "minus_treated",
    m = "primarily_leblango",
    y = "EL_EGRA_PCA_Index"))
})

