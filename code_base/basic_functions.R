#------------------------------------------------------------------------------------------------
#             Basic Functions
#------------------------------------------------------------------------------------------------


#' Sample from Symmetric Bernoulli Matrix
#' 
#' Samples a symmetric Bernoulli matrix with indendepent entries elementwise where 
#' the success probabilities are the entries of P.
#'
#' @param P the symmetric probability matrix containing success probabilities
#' 
#' @return A the sampled symmetric bernoulli matrix 
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
sampP <- function(P){
  #set up matrix
  n <- ncol(P)
  A <- matrix(0, n, n)
  
  #sample from A 
  A[upper.tri(A)] <- 1*(runif(sum(upper.tri(P))) < P[upper.tri(P)])
  A <- A + t(A)
  return(A)
}

#' @describeIn Make Omnibus Matrix
#' 
#' Creates the full nm x nm omnibus matrix
#'
#' @param mats a list of m, n x n matrices
#' 
#' @return Atilde the omnibus matrix of matrices in mats
#' @author Benjamin Draves <dravesb@bu.edu>
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

#' @describeIn Make Centered Omnibus Matrix
#' 
#' Creates the full nm x nm omnibus matrix where each n x n block is centered by the sample average in mats
#'
#' @param mats a list of m, n x n matrices
#' 
#' @return The omnibus matrix of matrices in mats centered by the sample mean of matrices in mats
#' @author Benjamin Draves <dravesb@bu.edu>
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
script <- RCurl::getURL("https://raw.githubusercontent.com/youngser/gmmase/master/R/getElbows.R", ssl.verifypeer = FALSE)
eval(parse(text = script))

#' @describeIn Function to perform graph adjacency spectral embedding (ASE)
#' 
#' @param A adjacency matrix
#' @param d number of joint embedding dimensions. If NA, dimension is chosen automatically
#' @param d.max maximum number of embedding dimensions to try when d is not provided. Default is sqrt(ncol(A)).
#' @param diag.augment whether to do diagonal augmentation (TRUE/FALSE)
#' @param elbow number of elbow selected in Zhu & Ghodsi method. Default is 1.
#' @return A matrix with n rows and d columns containing the estimated latent positions
#' 
#' @author Jes\'us Arroyo <jesus.arroyo@jhu.edu>
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

#' @describeIn Function to perform the indefinite spectral embedding (ISE) by choosing largest eigenvalues in magnitude.
#' 
#' @param A adjacency matrix
#' @param d number of embedding dimensions chosen from both sides of the spectral. If NA, dimension is chosen automatically
#' @param d.max maximum number of embedding dimensions to try when d is not provided. Default is sqrt(ncol(A)).
#' @param diag.augment whether to do diagonal augmentation (TRUE/FALSE)
#' @param elbow number of elbow selected in Zhu & Ghodsi method. Default is 1.
#' @return A matrix with n rows and d columns containing the estimated latent positions from the
#' 
#' @author Benjamin Draves <dravesb@bu.edu>
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

#' @describeIn Fetches graph-specific latent position estimates from Omnibus embedding
#' 
#' Subsets the nm x d omnibus embedding into a m - list of n x d 
#' matrices corresponding to graph specific latent position estimates
#'
#' @param Lhat the omnibus embedding matrix
#' @param m the numbe of graphs
#' 
#' @return a list of length m with each entry a n x d matrix corresponding to latent position estimates of each graph 
#' @author Benjamin Draves <dravesb@bu.edu>
get_Xhat_list <- function(Lhat, m){
  #fetch m 
  n <- nrow(Lhat) / m 
  
  #make X list
  Xhat_list <- list()
  for(g in 1:m) Xhat_list[[g]] <- Lhat[(n*(g-1) + 1):(n*g), ]
  
  #return list
  return(Xhat_list)
}
  
