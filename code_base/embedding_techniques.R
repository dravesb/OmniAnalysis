#---------------------------------------
#     Embedding Techniques
#
# omnibus and mase-tls
#---------------------------------------

#MASE and Omni functions
omni <- function(A_list, d, ...){
  #make A til 
  Atil <- make_omni(A_list)
  
  #embedd 
  Z_hat <- ase(Atil, d)
  
  #break into graph latent positions
  X_list <- list(); n <- nrow(A_list[[1]])
  for(g in 1:length(A_list)){
    start <- n*(g-1) + 1
    stop <- n*g
    X_list[[g]] <- as.matrix(Z_hat[start:stop,], nrow = n)
  }
  
  #return 
  return(X_list)
  
  
}
mase <- function(A_list, d_list, d, ...){
  #fetch individual graph ase
  X_ase <- lapply(1:length(A_list), function(g) ase(A_list[[g]], d_list[g]))
  
  #construct X
  X <- do.call('cbind', X_ase)
  
  #TLS
  X.svd <- svd(X, nu = d, nv = d)
  X_hat <- X.svd$u %*% tcrossprod(diag(c(X.svd$d[1:d]), nrow = d, ncol = d), X.svd$v)
  
  #break into list of latent positions
  X_list <- list()
  for(i in 1:length(d_list)){
    start <- cumsum(d_list)[i] - d_list[i] + 1
    stop <- cumsum(d_list)[i]
    X_list[[i]] <- as.matrix(X_hat[, start:stop], ncol = d_list[i])
  }
  
  #return 
  return(X_list)
}