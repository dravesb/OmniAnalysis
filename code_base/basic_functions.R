#------------------------------------------------------------------------------------------------
#             Basic Functions
#
# sampP: samples adjacency matrix from P 
# make_omni: takes list of adjacency matrix and makes Omnibus matrix
# make_centered_omni: takes list of adjacency matrix and makes Centered Omnibus matrix
# ase: Adjacency Spectral Embedding
# ise: Indefinite Spectral Embedding
# get_Xhat_list: Transform Lhat into list of X matrices
#------------------------------------------------------------------------------------------------

#sampling functions
sampP <- function(P) {
  #set up matrix
  n = ncol(P)
  A = matrix(0, n, n)
  
  #sample from A 
  A[upper.tri(A)] <- 1*(runif(sum(upper.tri(P))) < P[upper.tri(P)])
  A = A + t(A)
  return(A)
}

#Make Omni function
make_omni <- function(mats){
  #H(x) = (1x^T + x1^T)/2
  H <- function(g,m){
    ones <- rep(1, m)
    e <- diag(m)[,g]
    .5 * (tcrossprod(ones,e) + tcrossprod(e,ones))
  }
  
  #sum up each kronecker 
  m <- length(mats)
  Reduce("+", lapply(1:m, function(x) kronecker(H(x,m), mats[[x]])))
}

make_centered_omni <- function(mats){
  #H(x) = (1x^T + x1^T)/2 - 1/m J_{mm}
  H <- function(g,m){
    ones <- rep(1, m)
    e <- diag(m)[,g]
    .5 * (tcrossprod(ones,e) + tcrossprod(e,ones)) - 1/m * matrix(1, m, m)
  }
  
  #sum up each kronecker 
  m <- length(mats)
  Reduce("+", lapply(1:m, function(x) kronecker(H(x,m), mats[[x]])))
}

#source get Elbows function
script <- getURL("https://raw.githubusercontent.com/youngser/gmmase/master/R/getElbows.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#ASE Embedding
ase <- function(A, d = NA, d.max = ncol(A), diag.augment = T, elbow = 1) {
  require(rARPACK)
  # Diagonal augmentation
  if(diag.augment & sum(abs(diag(A))) == 0) {
    deg = colSums(A)
    n = ncol(A)
    diag(A) = deg / (n-1)
  }
  if(is.na(d)) {
    #Want to get the PD part of A
    eig <- eigs(as(A, "dgeMatrix"), d.max, which = 'LR')
    vals <- sort(x =  eig$values[which(eig$values > 0)], decreasing = TRUE) 
    d <- getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(eig$values >= vals[d])
    V <- eig$vectors[,selected.eigs, drop = F]
    D <- diag(sqrt(abs(eig$values[selected.eigs])), nrow = d)
    X <- V %*% D
    return(X)
  } else {
    eig <- eigs(as(A, "dgeMatrix"), k = d, which = 'LR')
    X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
    return(X) 
  }
}

#ISE Embedding
ise <- function(A, d = NA, d.max = ncol(A), diag.augment = T, elbow = 1) {
  require(rARPACK) 
  if(is.na(d)) {
    #augment diagonal with row means
    if(diag.augment & sum(abs(diag(A))) == 0) {
      deg = colSums(A)
      n = ncol(A)
      diag(A) = deg / (n-1)
    }
    #find embedding dimension and embedd
    eigv <- eigs(as(A, "dgeMatrix"), d.max)
    vals <- sort(x =  abs(eigv$values), decreasing = TRUE)
    d <- getElbows(vals, plot = F)[elbow]
    selected.eigs <- which(abs(eigv$values) >= vals[d])
    X <- eigv$vectors[,selected.eigs, drop = F] %*% diag(sqrt(abs(eigv$values[selected.eigs])), nrow = d)
    D <- sign(eigv$values[selected.eigs])
  } else{
    eig <- eigs(as(A, 'dgeMatrix'), k = d) #defaults to chooseing largest in magnitude eigenvalues
    X <- eig$vectors %*% diag(sqrt(abs(eig$values)), nrow = d)
    D <- sign(eig$values) 
  }
  return(list(X=X, D=D))
}

#Get Xhat matrices
get_Xhat_list <- function(Lhat, m){
  #fetch m 
  n <- nrow(Lhat) / m 
  
  #make X list
  Xhat_list <- list()
  for(g in 1:m) Xhat_list[[g]] <- Lhat[(n*(g-1) + 1):(n*g), ]
  
  #return list
  return(Xhat_list)
}
  
