library(ape)
library(tools)
library(stringr)

load_data <- function(data_file_name){
    T = load_tree(data_file_name)

    ## number of leaves in every gene tree
    n <- length(T[[1]]$tip.label) 
    # dimension of tree space
    e <- choose(n, 2)
    
    N <- length(T) ## sample size of the data set 
    D <- matrix(rep(0, N*e), N, e) 

    for(i in 1:N){
        u <- T[[i]]
        # u$edge.length <- u$edge.length/max(u$edge.length)
        dist_mat <- cophenetic(u)
        labs <- str_sort(u$tip.label)
        dist_mat <- dist_mat[labs,labs]
        dist_mat <- dist_mat / max(dist_mat)
        D[i, ] <- dist_mat[lower.tri(t(dist_mat))]
    }
    return(D)
}

load_tree <- function(data_file_name){
  file_type = file_ext(data_file_name)
  if(file_type %in% c("nex","dat") ){
    T = ape::read.nexus(data_file_name)
  } 
  if (file_type %in% c("tre","txt") ){
    T = ape::read.tree(data_file_name)
  }
  T = c(T)
  for(i in 1:length(T)){
    u <- T[[i]]
    # u$edge.length <- u$edge.length/max(u$edge.length)
    u$tip.label <- tolower(u$tip.label) 
    T[[i]] = u
  }
  T
}
