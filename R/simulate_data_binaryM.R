#Write a function for simulating data given the partial densities
simulate_data_binaryM <-
function(p_y_m1d1,
         p_y_m1d0, 
         p_y_m0d1, 
         p_y_m0d0, 
         p_m_1,
         p_m_0,
         p_d = 0.5,
         n, 
         seed = NULL, 
         yvalues = 1:length(p_y_m1d1)){
  
  #Set the seed if it is provided
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  #Draw d as a binomial with probability p_d
  d <- rbinom(n = n, size = 1, prob = p_d)

  #Depending on the value of d, draw m with probability p_m_1 or p_m_0
  m <- rbinom(n = n, size = 1, prob = d*p_m_1 + (1-d)*p_m_0)

  #Depending on the value of m and, draw y from yvalues with the probabilities given in p_y_m1d1, p_y_m1d0, p_y_m0d1, p_ym0d0
  y <- rep(NA, n)
  y[d ==1 & m==1] <- sample(x = yvalues, size = sum(d==1 & m==1), replace = TRUE, prob = p_y_m1d1)
  y[d ==1 & m==0] <- sample(x = yvalues, size = sum(d==1 & m==0), replace = TRUE, prob = p_y_m0d1)
  y[d ==0 & m==1] <- sample(x = yvalues, size = sum(d==0 & m==1), replace = TRUE, prob = p_y_m1d0)
  y[d ==0 & m==0] <- sample(x = yvalues, size = sum(d==0 & m==0), replace = TRUE, prob = p_y_m0d0)
 
  #Create a data frame
  df <- data.frame(d = d, m = m, y = y)
  
  return(df)
}
  
