rm(list=ls())
source("coalescent_model/load_data.R")

e=45

trop_dist <- function(x,y){
  max(x-y) - min(x-y)
}

eucl_dist <- function(x,y){
  d = x - y
  sqrt(c(d) %*% c(d))
}

dir = "coalescent_model/tree_sims/gen_error_data/"
Rs = c(.1,.2,.4,.6,.8,1,1.2,1.4,1.6,1.8,2,4,6,8,10)

S1 = load_data(paste(dir,"species_tree_1.nex",sep=""))
S2 = load_data(paste(dir,"species_tree_2.nex",sep=""))
species.dist = trop_dist(S1,S2)

library(distory)
bhv.dist <- function(s,u){
  x = c(s,u)
  dist.multiPhylo(x)[1]
}


sigma_est <- function(S1,u1,S2,u2){
  N = dim(u1)[1]
  distances.trop1 = rep(0,N)
  for ( i in 1:N ) distances.trop1[i] = trop_dist(S1,u1[i,])
  distances.trop2 = rep(0,N)
  for ( i in 1:N ) distances.trop2[i] = trop_dist(S2,u2[i,])
  hist(distances.trop1)
  
  distances.trop = c(distances.trop1,distances.trop2)

 ( mean(distances.trop1) + mean(distances.trop2)) /(2 * (e-1))

}

h <- function(x,dist) dist(S1,x) - dist(S2,x)

prob_est <-  function(S1,u1,S2,u2,dist){
  is.phylo = class(u1[[1]]) == "phylo" 
  N = ifelse(is.phylo,length(u1),  dim(u1)[1])
  h1 = rep(0,N)
  h2 = rep(0,N)
  if(is.phylo){
    for ( i in 1:N ) h1[i] = h(u1[[i]],dist)
    for ( i in 1:N ) h2[i] = h(u2[[i]],dist)
  } else {
    for ( i in 1:N ) h1[i] = h(u1[i,],dist)
    for ( i in 1:N ) h2[i] = h(u2[i,],dist)
  }
  (mean(h1>0) + mean(h2<0))/2
}

sigmas = list()
probs_trop = list()
probs_eucl = list()
for(R in Rs){
  print(R)
  u1 = load_data(paste(dir,"gene_trees_",R,"_1.nex",sep=""))
  u2 = load_data(paste(dir,"gene_trees_",R,"_2.nex",sep=""))
  sigmas[[paste(R)]] = sigma_est(S1,u1,S2,u2)
  probs_trop[[paste(R)]] = prob_est(S1,u1,S2,u2,trop_dist)
  probs_eucl[[paste(R)]] = prob_est(S1,u1,S2,u2,eucl_dist)
}
S1 = load_tree(paste(dir,"species_tree_1.nex",sep=""))
S2 = load_tree(paste(dir,"species_tree_2.nex",sep=""))
probs_bhv = list()
for(R in Rs){
  print(R)
  u1 = load_tree(paste(dir,"gene_trees_",R,"_1.nex",sep=""))
  u2 = load_tree(paste(dir,"gene_trees_",R,"_2.nex",sep=""))
  probs_bhv[[paste(R)]] = prob_est(S1,u1,S2,u2,bhv.dist)
}

# plot(Rs,sigmas)
# plot(Rs,probs_trop,log="x")
# lines(Rs,probs_eucl)

library(expint)
upper_bound <- function(sigma){
  x = species.dist/(2*sigma)
  gammainc(e-1,x)/(2*gamma(e-1))
}

theoretical_probs = sapply(sigmas,upper_bound)

# plot(sigmas,probs,ylim=c(0,0.5))
# lines(sigma_grid,theoretical_probs)

sigmas = unlist(sigmas)
probs_bhv  = unlist(probs_bhv)
probs_trop = unlist(probs_trop)
probs_eucl = unlist(probs_eucl)

save(list=c("Rs","probs_eucl","probs_trop","probs_bhv","sigmas","theoretical_probs"),file="gen_error_data.RData")
