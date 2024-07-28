compute_bootstrap_draws <- function(f, df, d, m, y, numdraws = 100){
  n <- NROW(df)

  bootstrap_one_seed <- function(seed){
    set.seed(seed)
    obs <- sample.int(n=NROW(df),size =NROW(df),replace = TRUE)
    df_bootstrap <- df[obs,]
    return(f(df_bootstrap, d, m, y))
  }

  bootstrapDraws <- purrr::map_dfr(.x = 1:numdraws, .f = bootstrap_one_seed)
  return(bootstrapDraws)
}


compute_bootstrap_draws_clustered <- function(f, df, d, m, y,
                                              cluster = NULL,
                                              numdraws = 100,
                                              return_df = T,
                                              fix_n1 = T){
  n <- NROW(df)

  if(is.null(cluster)){
    ncluster <- n
    uniqueClusters <- 1:n
    df$cluster <- 1:n
    cluster <- "cluster"
  }else{
    uniqueClusters <- unique(df[[cluster]])
    ncluster <- length(uniqueClusters)
    df$cluster <- df[[cluster]]
  }

  #Compute the number of treated and control clusters
  # If fix_n1, then this is treated as fixed and we bootstrap separately from each
  treated_clusters <- unique(df$cluster[df[[d]] == 1] )
  untreated_clusters <- unique(df$cluster[df[[d]] == 0] )
  ncluster_1 <- length(treated_clusters)
  ncluster_0 <- length(untreated_clusters)

  bootstrap_oneseed <- function(seed){
    set.seed(seed)

    if(!fix_n1){
      #If we don't fix n1, we just draw a bootstrap sample of clusters
      bs_clusters <- sample(x=uniqueClusters,size=ncluster,replace = TRUE)
    
    }
    else {
      #If we fix n1, we draw a bootstrap sample of clusters for treated/control, then combine
      bs_clusters1 <- sample(x=treated_clusters,size=ncluster_1,replace = TRUE)
      bs_clusters0 <- sample(x=untreated_clusters,size=ncluster_0,replace = TRUE)
      bs_clusters <- c(bs_clusters0, bs_clusters1)
    
    }
    
    bs_cluster_df <- data.frame(blah = 1:ncluster)
    bs_cluster_df[[cluster]] <- bs_clusters
    bs_cluster_df$blah <- NULL

    df_bs <- dplyr::left_join(bs_cluster_df,
                           df,
                           by = cluster)
    return(f(df = df_bs, d = d,m = m, y = y))
  }
  
  # First check that the treatment status is the same within each cluster
  if ((!is.null(cluster)) & (fix_n1)) {
    
    cluster_list <- unique(df[[cluster]])
    
    for (i in 1:length(cluster_list)) {
      temp <- df[ df[[cluster]] == cluster_list[i], d ]
      same_treat_status <- (sum(temp) == 0) | (sum(temp) == nrow(temp))
      
      if(!same_treat_status) {
        stop("You have set fix_n1 = TRUE, which fixes the number of treated clusters in the bootstrap. 
             To use this option, treatment status must be constant within clusters. 
             Use fix_n1 = FALSE if treatment status varies within clusters")
      }
    }
  }
  
  
  #If return_df, we return a df that binds that rows of the bootstrap draws
  #Otherwise, we return a list
  if(return_df){
  bootstrapDraws <- purrr::map_dfr(.x = 1:numdraws, .f = bootstrap_oneseed)
  }else{
  bootstrapDraws <- purrr::map(.x = 1:numdraws, .f = bootstrap_oneseed)

  }
  return(bootstrapDraws)
}

get_bootstrap_draw_clustered <- function(df, d, m, y, cluster, fix_n1 = F){
  n <- NROW(df)

  if(is.null(cluster)){
    ncluster <- n
    uniqueClusters <- 1:n
    df$cluster <- 1:n
    cluster <- "cluster"
  }else{
    uniqueClusters <- unique(df[[cluster]])
    ncluster <- length(uniqueClusters)
    df$cluster <- df[[cluster]]
  }

  # Compute the number of treated and control clusters
  # If fix_n1, then this is treated as fixed and we bootstrap separately from each
  treated_clusters <- unique(df$cluster[df[[d]] == 1] )
  untreated_clusters <- unique(df$cluster[df[[d]] == 0] )
  ncluster_1 <- length(treated_clusters)
  ncluster_0 <- length(untreated_clusters)

  if(!fix_n1){
    #If we don't fix n1, we just draw a bootstrap sample of clusters
    bs_clusters <- sample(x=uniqueClusters,size=ncluster,replace = TRUE)
  }else{
    #If we fix n1, we draw a bootstrap sample of clusters for treated/control, then combine
    bs_clusters1 <- sample(x=treated_clusters,size=ncluster_1,replace = TRUE)
    bs_clusters0 <- sample(x=untreated_clusters,size=ncluster_0,replace = TRUE)
    bs_clusters <- c(bs_clusters0, bs_clusters1)
  }
  bs_cluster_df <- data.frame(blah = 1:ncluster)
  bs_cluster_df[[cluster]] <- bs_clusters
  bs_cluster_df[[d]] <- 1
  bs_cluster_df[[d]][1:ncluster_0] <- 0
  bs_cluster_df$blah <- NULL

  df_bs <- dplyr::left_join(bs_cluster_df,
                            df[,c(m, y, cluster)],
                            by = cluster)

  return(df_bs[,c(d, m, y)])
}
