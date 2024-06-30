test_that("Reasonable arguments should run while nonsensical arguments raise error", {
  
  ### Use binary when M is not binary
  df <- burstzyn_data
  df$signed_up_number <- sample(0:3, nrow(df), replace = TRUE)
  expect_error(test_sharp_null(
    df = df,
    d = "condition2",
    m = "signed_up_number",
    y = "index", 
    method = "CS",
    use_binary = TRUE))

})

test_that("Expect null to be / not be rejected in specific scenarios", {
  
  methods <- c("CS", "ARP")
  
  ### Do not reject when treated and control outcomes are the same
  df <- burstzyn_data
  df_treated <- df[df$condition2 == 1, ]
  df_treated$condition2 <- 0
  df_test <- rbind(df[df$condition2 == 1, ], df_treated)
  
  for (method in methods) {
    result <- test_sharp_null(
      df = df_test,
      d = "condition2",
      m = "signed_up_number",
      y = "index", 
      method = method)
    expect_true(result$reject == FALSE)
  }
  
  ### Do not reject when Y does not depends on D, binary M
  N <- 300
  d <- sample(0:1, N, replace = TRUE)
  m <- sample(0:1, N, replace = TRUE)
  m[d == 1] <- 1
  y <- 2*m + rnorm(N, mean = 0, sd = 0.5)
  df_sim <- data.frame(treated = d,
                   mediator = m,
                   outcome = y)
  
  for (method in methods) {
    result <- test_sharp_null(
      df = df_sim,
      d = "treated",
      m = "mediator",
      y = "outcome", 
      method = method)
    expect_true(result$reject == FALSE)
  }
  
  ### Should reject when Y depends on both M and D, binary M
  N <- 300
  d <- sample(0:1, N, replace = TRUE)
  m <- sample(0:1, N, replace = TRUE)
  # m[d == 1] <- 1
  y <- 2*m + 30*d + rnorm(N, mean = 0, sd = 0.5)
  df_sim <- data.frame(treated = d,
                       mediator = m,
                       outcome = y)
  
  for (method in methods) {
    result <- test_sharp_null(
      df = df_sim,
      d = "treated",
      m = "mediator",
      y = "outcome", 
      method = method)
    expect_true(result$reject == TRUE)
  }
  
})
