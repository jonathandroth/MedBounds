test_that("Cox and Shi no nuisance function gives reasonable output", {

  #One moment binding, large
  cs <- cox_shi_nonuisance(Y = c(2,-1), sigma = diag(2))
  expect_equal(as.numeric(cs$test_stat), 2^2)
  expect_equal(as.numeric(cs$cv), qchisq(p=0.95, df = 1))
  expect_true(as.logical(cs$reject))

  #One moment binding, small
  cs <- cox_shi_nonuisance(Y = c(0.5,-1), sigma = diag(2))
  expect_equal(as.numeric(cs$test_stat), 0.5^2)
  expect_equal(as.numeric(cs$cv), qchisq(p=0.95, df = 1))
  expect_false(as.logical(cs$reject))

  #Two moments binding, large
  cs <- cox_shi_nonuisance(Y = c(5,6), sigma = diag(2))
  expect_equal(as.numeric(cs$test_stat), 5^2 + 6^2)
  expect_equal(as.numeric(cs$cv), qchisq(p=0.95, df = 2))
  expect_true(as.logical(cs$reject))

  #No moments binding
  cs <- cox_shi_nonuisance(Y = c(-5,-6), sigma = diag(2))
  expect_equal(as.numeric(cs$test_stat), 0)
  expect_false(as.logical(cs$reject))

  #Non-full rank variance
  #No moments binding
  cs <- cox_shi_nonuisance(Y = c(-1,-2), sigma = matrix(c(1,1,1,1),nrow = 2))
  expect_equal(as.numeric(cs$test_stat), 0)
  expect_false(as.logical(cs$reject))

  #Non-full rank variance
  #One moment violated
  cs <- cox_shi_nonuisance(Y = c(2,-1),
                           sigma = matrix(c(1,1,1,1),nrow = 2))
  expect_equal(as.numeric(cs$test_stat), 2^2)
  expect_equal(as.numeric(cs$cv), qchisq(p=0.95, df = 1))
  expect_true(as.logical(cs$reject))


})
