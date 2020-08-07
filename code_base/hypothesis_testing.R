#------------------------------------------------------------------------------------------------
#             Hypothesis Testing
#
# SigmaD_hat - estimation asympoptic variance for test statistic
# W_hat_stat -returns the two - graph test statistic and corresponding p-value
# T_hat_stat - returns the test statistic of Levin et al. 
# get_test_stat - return one of the two test statistics
#------------------------------------------------------------------------------------------------

#variance functions
Sigma_D_hat <- function(L_hat){
  #locally define helper function
  Sigma_sum_hat <- function(Xbar, vertex){
    #fetch n, d
    n <- nrow(Xbar)
    d <- ncol(Xbar)
    
    #get inner sum estimate
    total <- diag(0, nrow = d, ncol = d)
    for(j in (1:n)[-vertex]){
      total <- total + as.numeric(crossprod(Xbar[vertex,], Xbar[j,]) - crossprod(Xbar[vertex,], Xbar[j,])^2) * tcrossprod(Xbar[j,])
    }
    
    #return sum 
    return(1/(n-1) * total)
    
  }
  
  #fetch n,m,d
  d <- ncol(L_hat)
  m <- 2
  n <- nrow(L_hat)/m
  
  #set S2_Delta_hat and Omnibar
  S2_Delta_hat <- crossprod(L_hat)/(n*m)
  Omnibar <- 0.5*(L_hat[1:n,] + L_hat[(n+1):(2*n),])
  
  #set up variance list 
  var.here <- list()
  for(i in 1:n){
    var.here[[i]] <- 0.5 * solve(S2_Delta_hat) %*% Sigma_sum_hat(Omnibar,i) %*% solve(S2_Delta_hat)
  }
  
  #return estimated variances for each row of D
  return(var.here)
  
  
}

#W_hat test statistic test statistics
W_hat_stat <- function(L_hat, return_p_val = FALSE, correction = 0){
  #locally define helper function
  psd_proj <- function(M){
    #decomp 
    eigen_decomp <- eigen(M)
    U <- eigen_decomp$vectors
    S <- eigen_decomp$values
    
    #threshold eigenvalues
    S.thres <- diag(sapply(S, function(x) max(c(0,x))))
    
    #return
    return(tcrossprod(U %*% S.thres, U))
    
  }
  #fetch n, m ,d
  m <- 2
  n <- nrow(L_hat)/m
  d <- ncol(L_hat)
  
  #get D_hat
  D_hat <- L_hat[1:n, ] - L_hat[(n+1):(2*n), ]
  
  #fetch vertex level variance
  var.here <- Sigma_D_hat(L_hat)
  
  #project onto psd cone and then psuedo inverse
  var.here.invert <- lapply(var.here, function(X) ginv(psd_proj(X)))
  
  #get vertex level test stat
  Wi_hat_list <- lapply(1:n, function(i) n*crossprod(D_hat[i,], var.here.invert[[i]]) %*% D_hat[i,])
  
  #get W_hat
  W_hat <- as.numeric(Reduce('+', Wi_hat_list))
  
  #return vector
  if(return_p_val){
    p_val <- pchisq(W_hat, n*d + correction, lower.tail = FALSE)
    ret <- c(W_hat, p_val); names(ret) <- c('Test Stat', 'p-value')
    return(ret)  
  }else{
    return(W_hat)
  }
  
  
}

#T_hat test statistic of Levin et al. 
T_hat_stat <- function(L_hat){
  #fetch n
  m <- 2
  n <- nrow(L_hat)/m
  
  #get D_hat
  D_hat <- L_hat[1:n, ] - L_hat[(n+1):(2*n), ]
  
  #return statistic 
  T_hat <- norm(D_hat, type = 'F')
  return(T_hat)
  
}

#get statistic: option 1 for covariance corrected 
#               option 2 for Frobenius distance
get_test_stat <- function(L_hat, option = 1, ...) return(ifelse(option == 1, W_hat_stat(L_hat, ...), T_hat_stat(L_hat)))
