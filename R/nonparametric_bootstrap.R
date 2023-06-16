compute_bootstrap_draws <- function(f, df, d, m, y, w = NULL, numdraws = 100){
  n <- NROW(df)

  bootstrap_one_seed <- function(seed){
    set.seed(seed)
    obs <- sample.int(n=NROW(df),size =NROW(df),replace = TRUE)
    df_bootstrap <- df[obs,]
    return(f(df_bootstrap, d, m, y, w = "dirichletDraws"))
  }

  bootstrapDraws <- purrr::map_dfr(.x = 1:numdraws, .f = bootstrap_one_seed)
  return(bootstrapDraws)
}


compute_bootstrap_draws_clustered <- function(f, df, d, m, y, cluster = NULL, numdraws = 100){
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



  bootstrap_oneseed <- function(seed){
    set.seed(seed)
    bs_clusters <- sample(x=uniqueClusters,size=ncluster,replace = TRUE)

    bs_cluster_df <- data.frame(blah = 1:ncluster)
    bs_cluster_df[[cluster]] <- bs_clusters
    bs_cluster_df$blah <- NULL

    df_bs <- dplyr::left_join(bs_cluster_df,
                           df,
                           by = cluster)
    return(f(df_bs, d, m, y))
  }

  bootstrapDraws <- purrr::map_dfr(.x = 1:numdraws, .f = bootstrap_oneseed)
  return(bootstrapDraws)
}
