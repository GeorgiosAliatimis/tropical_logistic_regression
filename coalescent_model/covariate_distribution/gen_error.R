rm(list=ls())
source("coalescent_model/load_data.R")

e=45

trop_dist <- function(x,y){
  max(x-y) - min(x-y)
}

dir = "coalescent_model/tree_sims/gen_error_data/"
Rs = c(.1,.2,.4,.6,.8,1,1.2,1.4,1.6,1.8,2,4,6,8,10)

S1 = load_data(paste(dir,"species_tree_1.nex",sep=""))
S2 = load_data(paste(dir,"species_tree_2.nex",sep=""))
species.dist = trop_dist(S1,S2)

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

h <- function(x) trop_dist(S1,x) - trop_dist(S2,x)

prob_est <-  function(S1,u1,S2,u2){
  N = dim(u1)[1]
  h1 = rep(0,N)
  for ( i in 1:N ) h1[i] = h(u1[i,])
  h2 = rep(0,N)
  for ( i in 1:N ) h2[i] = h(u2[i,])
  (mean(h1>0) + mean(h2<0))/2
}

sigmas = list()
probs = list()
for(R in Rs){
  u1 = load_data(paste(dir,"gene_trees_",R,"_1.nex",sep=""))
  u2 = load_data(paste(dir,"gene_trees_",R,"_2.nex",sep=""))
  sigmas[[paste(R)]] = sigma_est(S1,u1,S2,u2)
  probs[[paste(R)]] = prob_est(S1,u1,S2,u2)
}
plot(Rs,sigmas)
plot(Rs,probs)


library(expint)
upper_bound <- function(sigma){
  x = species.dist/(2*sigma)
  gammainc(e-1,x)/(2*gamma(e-1))
}

sigma_grid = seq(0,0.05,length=1001)
theoretical_probs = sapply(sigma_grid,upper_bound)

plot(sigmas,probs,ylim=c(0,0.5))
lines(sigma_grid,theoretical_probs)
