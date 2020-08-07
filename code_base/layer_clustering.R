#------------------------------------------------------------------------------------------------
#           Layer Clustering and Visualizations
#
# dissim_mat - a function that constructs a dissimilarity matrix between layers 
# cluster_mat - preforms clustering on the dissimilarity matrix
#------------------------------------------------------------------------------------------------

#Dissimalarity matrix construction
dissim_mat <- function(A_list, plot_matrix = FALSE, option = 1){
  #make omnibus matrix
  Atil <- make_omni(A_list)
  m <- length(A_list)
  n <- nrow(A_list[[1]])
  
  #embedd
  L_hat <- ase(Atil)
  
  #set up D matrix
  D <- matrix(0, ncol = m, nrow = m)
  for(i in 1:(m-1)){
    for(j in (i+1):m){
      #get indices 
      start_i <- n*(i-1) + 1; end_i <- n*i
      start_j <- n*(j-1) + 1; end_j <- n*j
      
      #store test statistic
      D[i,j] <- get_test_stat(L_hat[c(start_i:end_i, start_j:end_j), ], option) 
    }
  }
  
  #symmetrize & plot if necessary
  D <- D + t(D)
  if(plot_matrix == TRUE){
    dev.new()
    image(Matrix(D), col.regions = adjustcolor('red', alpha.f = 0.6))
  }
  #return
  return(D)
  
}

#clustering of layers via clustering the Dissimilarity matrix
cluster_layers <- function(D, plot = F, k = 2, option = 1){
  if(option == 1){
    hc <- hclust(as.dist(D))
    
    #want to do this in an unsupervised fashion
    groups <- cutree(hc, k)
    return(groups)
  }
}

