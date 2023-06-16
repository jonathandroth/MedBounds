compute_posterior_draws <- function(f, df, d, m, y, w = NULL, numdraws = 100){
  n <- NROW(df)

  compute_posterior_helper <- function(seed){
    set.seed(seed)
    gammaDraws <- stats::rgamma(n = n, shape = 1) #draw indep gammas
    dirichletDraws <- gammaDraws / sum(gammaDraws) #normalize to get Dirichlet
    df$dirichletDraws <- dirichletDraws
    return(f(df, d, m, y, w = "dirichletDraws"))
  }

  posteriorDraws <- purrr::map_dfr(.x = 1:numdraws, .f = compute_posterior_helper)
  return(posteriorDraws)
}


compute_posterior_draws_clustered <- function(f, df, d, m, y, cluster = NULL, w = NULL, numdraws = 100){
  n <- NROW(df)

  if(is.null(cluster)){
    ncluster <- n
    uniqueClusters <- 1:n
    df$cluster <- 1:n
  }else{
    uniqueClusters <- unique(df[[cluster]])
    ncluster <- length(uniqueClusters)
    df$cluster <- df[[cluster]]
  }



  compute_posterior_helper <- function(seed){
    set.seed(seed)
    gammaDraws <- stats::rgamma(n = ncluster, shape = 1) #draw indep gammas for each cluster
    dirichletDraws <- gammaDraws / sum(gammaDraws) #normalize to get Dirichlet

    df <- dplyr::left_join(df,
                            data.frame(cluster = uniqueClusters, dirichletDraws = dirichletDraws),
                            by = "cluster")
    return(f(df, d, m, y, w = "dirichletDraws"))
  }

  posteriorDraws <- purrr::map_dfr(.x = 1:numdraws, .f = compute_posterior_helper)
  return(posteriorDraws)
}
