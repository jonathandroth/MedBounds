test_that("simple discretization of y works", {
  expect_equal(as.numeric(discretize_y(c(1,2,3,4,5), numBins = 2)),
                  c(1,1,1,2,2))

  expect_equal(as.numeric(discretize_y(c(1,2,3,4,5,6), numBins = 3)),
               c(1,1,2,2,3,3))

  expect_equal(as.numeric(discretize_y(c(rep(1,5), seq(from=6,to=10)), numBins = 3)),
               c(rep(1,5), rep(2,2), rep(3,3) ))

})
