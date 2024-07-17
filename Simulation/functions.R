library(lpSolve)
generate_M_list_lam <- function(d,r,lambda_list, kappa=1.5){
  M_list <- vector("list", length = length(lambda_list))
  for (i in 1:length(lambda_list)){
    U <- svd(matrix(rnorm(d*r), ncol = r))$u
    V <- svd(matrix(rnorm(d*r), ncol = r))$u
    Sigma <- diag(rev(seq(lambda_list[i], kappa*lambda_list[i], length.out = r)), nrow = r)
    M_list[[i]] <- U%*%Sigma%*%t(V)
  }
  return(M_list)
}
generate_M_list <- function(d,K,r,lambda,kappa=1.2){
  M_list <- vector(mode = "list", K)
  for (i in 1:K){
    U <- svd(matrix(rnorm(d*r), ncol = r))$u
    V <- svd(matrix(rnorm(d*r), ncol = r))$u
    Sigma <- diag(rev(seq(lambda, kappa*lambda, length.out = r)), nrow = r)
    M_list[[i]] <- U%*%Sigma%*%t(V)
  }
  Delta = min(dist(t(sapply(M_list, as.vector))))
  return(list("M_list" = M_list,"Delta"= Delta))
}
generate_two_M <- function(d,r,lambda,Delta,kappa=1.2){
  M_list <- vector(mode = "list", 2)
  U <- svd(matrix(rnorm(d*r), ncol = r))$u
  V <- svd(matrix(rnorm(d*r), ncol = r))$u
  Sigma <- diag(rev(seq(lambda, kappa*lambda, length.out = r)), nrow = r)
  M_list[[1]] <- U%*%Sigma%*%t(V)
  M_list[[2]] <- U%*%diag(diag(Sigma) + sign(Delta)*sqrt(abs(Delta)/r))%*%t(V)
  return(M_list)
}
generate_two_M_o <- function(d,r,lambda,kappa=1.2){
  M_list <- vector(mode = "list", 2)
  U <- svd(matrix(runif(d*(2*r)), ncol = 2*r))$u
  V <- svd(matrix(runif(d*(2*r)), ncol = 2*r))$u
  Sigma <- diag(rev(seq(lambda, kappa*lambda, length.out = r)), nrow = r)
  M_list[[1]] <- U[,1:r]%*%Sigma%*%t(V[,1:r])
  M_list[[2]] <- U[,(r+1):(2*r)]%*%Sigma%*%t(V[,(r+1):(2*r)])
  Delta <- norm(M_list[[1]]-M_list[[2]],"F")
  return(list("M_list" = M_list,"Delta"= Delta))
}
generate_M_list_relaxed <- function(d,K,r,lambda,kappa=1.2){
  M_list <- vector(mode = "list", K)
  for (i in 1:(K-1)){
    U <- svd(matrix(rnorm(d*r), ncol = r))$u
    V <- svd(matrix(rnorm(d*r), ncol = r))$u
    Sigma <- diag(rev(seq(lambda, kappa*lambda, length.out = r)), nrow = r)
    M_list[[i]] <- U%*%Sigma%*%t(V)
  }
  U <- svd(matrix(rnorm(d*r), ncol = r))$u
  V <- svd(matrix(rnorm(d*r), ncol = r))$u
  Sigma <- diag(rev(seq(0.3, kappa*0.3, length.out = r)), nrow = r)
  M_list[[K]] <- U%*%Sigma%*%t(V)
  Delta = min(dist(t(sapply(M_list, as.vector))))
  return(list("M_list" = M_list,"Delta"= Delta))
}
generate_data <- function(M_list, label, sigma=1){
  K <- length(M_list)
  d <- dim(M_list[[1]])[1]
  n <- length(label)
  X_tsnr <- array(NA, dim = c(d,d,n))
  for (i in 1:n){
    X_tsnr[,,i] <- M_list[[label[i]]] + matrix(rnorm(d*d, sd = sigma), ncol = d)
  }
  data_tsr <- as.tensor(X_tsnr)
  return(data_tsr)
}
generate_SBM<-function(P_mat) 
{
  d <- dim(P_mat)[1]
  A <- matrix(NA,ncol=d,nrow=d)
  diag(A) <- 0
  for (i in 1:(d-1))
  {
    for (j in (i+1):d)
    {
      A[i,j]=rbinom(1,1,P_mat[i,j])
      A[j,i]=A[i,j]
    }  
  }
  return(A)
}
generate_P_list <- function(d, r, K, maxp, gap_ratio){
  P_list <- vector(mode = "list", K)
  p_vec <- seq(maxp*gap_ratio, maxp, length.out = K)
  q_vec <- 0.5*p_vec
  membership_mat <- matrix(NA, nrow = K, ncol = d)
  for (k in 1:K){
    membership_mat[k,] = sample(r, d, replace = TRUE)
  }
  for (k in 1:K){
    B=matrix(q_vec[k],r,r)
    diag(B)=p_vec[k]
    Z=matrix(0,nrow=d,ncol=r)
    for(i in 1:d){
      Z[i,membership_mat[k,i]]=1
    }
    P_list[[k]] <- Z%*%B%*%t(Z)
  }
  Delta = min(dist(t(sapply(P_list, as.vector))))
  return(list("P_list" = P_list,"Delta"= Delta))
}
generate_network_data <- function(P_list, label){
  K <- length(P_list)
  d <- dim(P_list[[1]])[1]
  n <- length(label)
  X_tsnr <- array(NA, dim = c(d,d,n))
  for (i in 1:n){
    X_tsnr[,,i] <- generate_SBM(P_list[[label[i]]])
  }
  data_tsr <- as.tensor(X_tsnr)
  return(data_tsr)
}
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
fnorm <- function(X){
  ### Tensor Frobenius Norm \|X\|_F ###
  return(sqrt(sum(X^2)))
}

lrlloyds <- function(data_tsr, s_init, r_list, relaxed = FALSE){
  K <- length(unique(s_init))
  n <- length(s_init)
  Mhat_list <- vector(mode = "list", length = K)
  shat <- s_init
  tmax <- ceiling(3*log(n))
  ## Iterations
  if (relaxed==TRUE){
    for (t in 1:tmax){
      for (k in 1:K) {
        mean_mat <- apply(data_tsr@data[,,shat==k], c(1,2), mean)
        mean_svd <- svd(mean_mat, nu = r_list[k], nv = r_list[k])
        if (r_list[k]==1){
          Mhat_list[[k]] <- mean_svd$d[1]* mean_svd$u %*% t(mean_svd$v)
        } else{
          Mhat_list[[k]] <- mean_svd$u %*% diag(mean_svd$d[1:r_list[k]]) %*% t(mean_svd$v)
        }
      }
      Mhat_list[[which.min(sapply(Mhat_list, norm, "2"))]] <- 0
      for (i in 1:n) {
        shat[i] <- which.min(sapply(Mhat_list, function(x){norm(data_tsr@data[,,i]-x, "F")}))
      }
    }
  }else{
    for (t in 1:tmax){
      for (k in 1:K) {
        mean_mat <- apply(data_tsr@data[,,shat==k], c(1,2), mean)
        mean_svd <- svd(mean_mat, nu = r_list[k], nv = r_list[k])
        if (r_list[k]==1){
          Mhat_list[[k]] <- mean_svd$d[1]* mean_svd$u %*% t(mean_svd$v)
        } else{
          Mhat_list[[k]] <- mean_svd$u %*% diag(mean_svd$d[1:r_list[k]]) %*% t(mean_svd$v)
        }
      }
      for (i in 1:n) {
        shat[i] <- which.min(sapply(Mhat_list, function(x){norm(data_tsr@data[,,i]-x, "F")}))
      }
    }
  }
  ## Calculate WSS
  Mhat_tsr <- array(NA, dim = data_tsr@modes)
  for (i in 1:n) {
    Mhat_tsr[,,i] <- Mhat_list[[shat[i]]]
  }
  size <- prod(data_tsr@modes)
  mydf <- 0
  for (k in 1:K){
    mydf <- mydf + r_list[k]*(data_tsr@modes[1]+data_tsr@modes[2]+1)
  }
  tot.wss <- sum((Mhat_tsr-data_tsr@data)^2)
  bic <- log(tot.wss/size)+log(size) / (size) * mydf
  return(list("shat"=shat, "tot.wss"=tot.wss,  "bic"=bic, "Mhat_list"=Mhat_list))
}
random_lrlloyds <- function(data_tsr, rank, true_label, nstart, relaxed){
  res_mat <- matrix(NA, nrow = nstart, ncol = 3)
  for (i in 1:nstart){
    s_init <- sample(c(1,2), n, replace = TRUE)
    res <- lrlloyds(data_tsr, s_init, rank, relaxed)
    res_mat[i,1] <- res$bic
    res_mat[i,2] <- effective_hamming_error(s_init,true_label)
    res_mat[i,3] <- effective_hamming_error(res$shat,true_label)
  }
  res_mat
}
lr_spec_init <- function(tsr, rU, rV, K){
  HOSVD_init <- hosvd(tsr, ranks = c(rU, rV, K))
  U_int <- HOSVD_init$U[[1]]
  V_int <- HOSVD_init$U[[2]]
  Ghat <- ttl(tsr, list(U_int%*%t(U_int),V_int%*%t(V_int)), c(1,2))
  M3Ghat <- k_unfold(Ghat, m=3)@data
  clustering_res <- kmeans(M3Ghat, K, nstart = 30)$cluster
  return(clustering_res)
}
spec_init <- function(tsr,K){
  Y <- t(k_unfold(tsr, m=3)@data)
  Uk <- svd(Y,nu=K)$u
  Yhat <- (Uk)%*%(t(Uk)%*%Y)
  clustering_res <- kmeans(t(Yhat), K, nstart = 20)
  return(list("est_label"=clustering_res$cluster,"center"=clustering_res$centers))
}
vec_lloyd <- function(tsr, K){
  n <- tsr@modes[3]
  spec_res <- spec_init(tsr,K)
  s_vec_lloyds <- kmeans(k_unfold(tsr, m=3)@data, spec_res$center, iter.max = 3*log(n),
                         algorithm = "Lloyd", nstart = 30)$cluster
  s_vec_init <- spec_res$est_label
  return(list("init"=s_vec_init, "lloyd"=s_vec_lloyds))
}
lrlloyds_iter_err <- function(data_tsr, s_init, r_list, true_label){
  K <- length(unique(s_init))
  n <- length(s_init)
  Mhat_list <- vector(mode = "list", length = K)
  shat <- s_init
  tmax <- ceiling(3*log(n))
  s_err <- rep(NA, tmax+1)
  s_err[1] <- effective_hamming_error(shat,true_label)
  cat(s_err[1],"\n")
  ## Iterations
  for (t in 1:tmax){
    for (k in 1:K) {
      mean_mat <- apply(data_tsr@data[,,shat==k], c(1,2), mean)
      mean_svd <- svd(mean_mat, nu = r_list[k], nv = r_list[k])
      if (r_list[k]==1){
        Mhat_list[[k]] <- mean_svd$d[1]* mean_svd$u %*% t(mean_svd$v)
      } else{
        Mhat_list[[k]] <- mean_svd$u %*% diag(mean_svd$d[1:r_list[k]]) %*% t(mean_svd$v)
      }
    }
    for (i in 1:n) {
      shat[i] <- which.min(sapply(Mhat_list, function(x){norm(data_tsr@data[,,i]-x, "F")}))
    }
    s_err[t+1] <- effective_hamming_error(shat,true_label)
    cat(s_err[t+1],"\n")
  }
  ## Calculate WSS
  return(s_err)
}
generate_data_X <- function(M,s,n,d){
  X_array <- array(0,dim=c(d,d,n))
  for (i in 1:n){
    Z <- matrix(rnorm(d*d), ncol = d)
    X_array[,,i] <- s[i]*M+Z
  }
  return(X_array)
}

spec_aggregation <- function(M, X_array, n,d,r){
  X_list <- vector("list", length = n)
  for (i in 1:n){
    X_list[[i]] <- X_array[,,i] 
  }
  Xf <- do.call(cbind,X_list)
  ## Spectral initialization
  u_1hat <- svd(Xf, nu=1)$u
  X_array_tu1 <- apply(X_array, 3, function(x){t(x)%*%u_1hat})
  v_1hat <- svd(X_array_tu1, nu=1)$u
  ## Spectral refinement
  u1Xv1 <- apply(X_array, 3, function(x){t(u_1hat)%*%x%*%v_1hat})
  X_matrix <- matrix(X_array, nrow = d*d)
  M_1 <- matrix(rowMeans(X_matrix%*%diag(u1Xv1)), ncol = d)-u_1hat%*%t(v_1hat)
  svd_M_1 <- svd(M_1, nu=r,nv=r)
  U_tilde <- svd_M_1$u
  V_tilde <- svd_M_1$v
  ## Aggregation
  trUXV <- apply(X_array, 3, function(x){sum(diag(t(U_tilde)%*%x%*%V_tilde))})
  M_r <- matrix(rowMeans(X_matrix%*%diag(trUXV)), ncol = d)-U_tilde%*%t(V_tilde)
  svd_M_r <- svd(M_r, nu=r,nv=r)
  M_check <- svd_M_r$u %*% diag(svd_M_r$d[1:r], nrow = r) %*%t(svd_M_r$v)
  Lambda_hat <- sqrt(max(mean(trUXV^2)-r,0))
  if (Lambda_hat==0){
    M_hat <- M_check*0; index <- 1
  } else{
    M_hat <- M_check/Lambda_hat; index <- 0
  }
  return(c(min(sqrt(sum((M_hat-M)^2)),sqrt(sum((M_hat+M)^2))), sqrt(sum(M^2)), Lambda_hat))
}
