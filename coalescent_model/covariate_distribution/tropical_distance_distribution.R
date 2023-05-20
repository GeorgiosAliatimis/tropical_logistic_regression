source("coalescent_model/load_data.R")
dir = "coalescent_model/data/Depth0.7/"
library("distory")
library("phangorn")

trop.dist <- function(x,y){
  max(x-y) - min(x-y)
}

eucl.dist <- function(x,y){
  sqrt( c(x-y) %*% c(x-y) )
}

bhv.dist <- function(s,u){
  x = c(s,u)
  dist.multiPhylo(x)[1]
}

S1 = load_data(paste(dir,"Species1.tre",sep=""))[[1]]
u1 = load_data(paste(dir,"Genes1_Depth0.7.nex",sep=""))
S2 = load_data(paste(dir,"Species2.tre",sep=""))[[1]]
u2 = load_data(paste(dir,"Genes2_Depth0.7.nex",sep=""))

for(dist_type in c("eucl","trop")){
  fun = eval(parse(text=paste(dist_type,".dist",sep="")))
  d1 = apply(u1,1,function(x) fun(x,S1))
  d2 = apply(u2,1,function(x) fun(x,S2))
  d  = c(d1,d2)
  write.csv(d,paste("./coalescent_model/covariate_distribution/",dist_type,"_dists.csv",sep=""))
  hist(d,prob=T,main=dist_type)
}

S1 = load_tree(paste(dir,"Species1.tre",sep=""))[[1]]
u1 = load_tree(paste(dir,"Genes1_Depth0.7.nex",sep=""))
S2 = load_tree(paste(dir,"Species2.tre",sep=""))[[1]]
u2 = load_tree(paste(dir,"Genes2_Depth0.7.nex",sep=""))

library(purrr)
for(dist_type in c("bhv","RF")){
  fun = eval(parse(text=paste(dist_type,".dist",sep="")))
  d1 = unlist(map(u1,function(x) fun(x,S1)))
  d2 = unlist(map(u2,function(x) fun(x,S2)))
  d  = c(d1,d2)
  write.csv(d,paste("./coalescent_model/covariate_distribution/",dist_type,"_dists.csv",sep=""))
  hist(d,prob=T,main=dist_type)
}