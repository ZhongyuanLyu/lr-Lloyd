library(lpSolve)
effective_hamming_error <- function(cluster1, cluster2) {
  mapping_to_list <- function(K, clust) {
    map <- vector("list", K)
    map[[1]] <- c(map[[1]], 1)
    ind <- 1
    for (i in 2:n) {
      for (j in 1:(i - 1)) {
        find <- 0
        if (clust[i] == clust[j]) {
          for (k in 1:ind) {
            if (j %in% map[[k]]) {
              map[[k]] <- c(map[[k]], i)
              find <- 1
              break
            }
          }
        }
        if (find == 1) break
        if (j == (i - 1)) {
          ind <- ind + 1
          map[[ind]] <- c(map[[ind]], i)
        }
      }
    }
    return(map)
  }
  clust1 <- unclass(as.ordered(cluster1))
  clust2 <- unclass(as.ordered(cluster2))
  if ((n <- length(clust1)) != length(clust2)) {
    warning("error: length not equal")
    return
  }
  if ((K <- length(table(clust1))) != length(table(clust2))) {
    warning("the number of clusters are not equal")
    return
  }
  list1 <- mapping_to_list(K, clust1)
  list2 <- mapping_to_list(K, clust2)
  cost_mat <- matrix(0, nrow = K, ncol = K)
  for (i in 1:K) {
    for (j in 1:K) {
      cost_mat[i, j] <- sum(!(list1[[i]] %in% list2[[j]]))
    }
  }
  error <- lp.assign(cost_mat)
  return(error$objval / n)
}
lloyds <- function(data_tsr, s_init){
  K <- length(unique(s_init))
  n <- length(s_init)
  Mhat_list <- vector(mode = "list", length = K)
  shat <- s_init
  tmax <- ceiling(3*log(n))
  for (t in 1:tmax){
    for (k in 1:K) {
      Mhat_list[[k]] <- apply(data_tsr@data[,,shat==k], c(1,2), mean)
    }
    for (i in 1:n) {
      shat[i] <- which.min(sapply(Mhat_list, function(x){norm(data_tsr@data[,,i]-x, "F")}))
    }
  }
  shat
}
lrlloyds <- function(data_tsr, s_init, r_list){
  K <- length(unique(s_init))
  n <- length(s_init)
  Mhat_list <- vector(mode = "list", length = K)
  shat <- s_init
  tmax <- ceiling(3*log(n))
  for (t in 1:tmax){
    for (k in 1:K) {
      mean_mat <- apply(data_tsr@data[,,shat==k], c(1,2), mean)
      mean_svd <- svd(mean_mat, nu = r_list[k], nv = r_list[k])
      if (r_list[k]==1){
        Mhat_list[[k]] <-mean_svd$d[1]* mean_svd$u %*% t(mean_svd$v)
      } else{
        Mhat_list[[k]] <- mean_svd$u %*% diag(mean_svd$d[1:r_list[k]]) %*% t(mean_svd$v)
      }
    }
    for (i in 1:n) {
      shat[i] <- which.min(sapply(Mhat_list, function(x){norm(data_tsr@data[,,i]-x, "F")}))
    }
  }
  shat
}