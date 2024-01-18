library(ape)
library(tools)
library(stringr)
library(MASS)

load_data <- function(T){
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
    labs <- str_sort(u$tip.label)[1:n]
    dist_mat <- dist_mat[labs,labs]
    dist_mat <- dist_mat / max(dist_mat)
    D[i, ] <- dist_mat[lower.tri(t(dist_mat))]
  }
  return(D)
}

# args = commandArgs(trailingOnly=TRUE)
# samplefreq=as.numeric(args[1])
# diagnfreq=as.numeric(args[2])
# trees_per_generation=diagnfreq/samplefreq
data_from_file <- function(file,ind){
    trees=read.nexus(file)
    # L = length(trees)%/%trees_per_generation
    # for(i in 1:L) {
    	# D=load_data(trees[(trees_per_generation*(i-1)+1):(trees_per_generation*i)])
      # print(dim(D))
    	# write.matrix(D,file=paste("v", i,".csv",sep=""))
    # }
    D = load_data(trees)
    print(dim(D))
    write.matrix(D,file=paste("v", ind,".csv",sep=""))
}

data_from_file("t1",1)
data_from_file("t2",2)
