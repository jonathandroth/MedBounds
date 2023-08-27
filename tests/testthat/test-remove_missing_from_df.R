test_that("don't remove anything if df is not missing", {
  df <- data.frame(y1 = c(1,1), m1 = c(1,0), d1 = c(0,1))
  df_processed <- remove_missing_from_df(df, m = "m1", d= "d1", y = "y1")
  expect_equal(df,df_processed)
})


test_that("remove missing rows from df property", {
  df <- data.frame(y1 = c(NA,1,1,1), m1 = c(1,NA,1,0), d1 = c(1,0,NA,1))
  df_processed <- remove_missing_from_df(df, m = "m1", d= "d1", y = "y1")
  expect_equal(df[4,],df_processed)
})


test_that("remove missing rows from df property", {
  df <- data.frame(y1 = c(1,1,1,1), m1 = c(1,NA,1,0),
                   d1 = c(1,0,NA,1),
                   w1 = c(NA,1,1,1))
  df_processed <- remove_missing_from_df(df, m = "m1", d= "d1", y = "y1", w = "w1")
  expect_equal(df[4,],df_processed)
})

