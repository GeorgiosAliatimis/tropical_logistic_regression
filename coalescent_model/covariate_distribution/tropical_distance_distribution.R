source("coalescent_model/load_data.R")
dir = "coalescent_model/data/Depth0.7/"

S1 = load_data(paste(dir,"Species1.tre",sep=""))
u1 = load_data(paste(dir,"Genes1_Depth0.7.nex",sep=""))
S2 = load_data(paste(dir,"Species2.tre",sep=""))
u2 = load_data(paste(dir,"Genes2_Depth0.7.nex",sep=""))

trop_dist <- function(x,y){
  max(x-y) - min(x-y)
}

eucl_dist <- function(x,y){
  sqrt( c(x-y) %*% c(x-y) )
}

N = dim(u1)[1]
e = dim(u1)[2]
distances.trop1 = rep(0,N)
for ( i in 1:N ) distances.trop1[i] = trop_dist(S1,u1[i,])
distances.trop2 = rep(0,N)
for ( i in 1:N ) distances.trop2[i] = trop_dist(S2,u2[i,])
hist(distances.trop1, prob=TRUE)
hist(distances.trop2, prob=TRUE)

distances.trop = c(distances.trop1,distances.trop2)


sigma_tr = ( mean(distances.trop1) + mean(distances.trop2)) /(2 * (e-1))
sigma_tr1 = mean(distances.trop1) / (e-1)
sigma_tr2 = mean(distances.trop2) / (e-1) 

l_trop = - 2*N*log( factorial(e) ) - (e-1) * N * (log(sigma_tr1)  + log(sigma_tr2)) - (e-1) * 2 * N 

# distances.eucl = rep(0,N)
# for ( i in 1:N ) distances.eucl[i] = eucl_dist(S,u[i,])

distances.eucl1 = rep(0,N)
for ( i in 1:N ) distances.eucl1[i] = eucl_dist(S1,u1[i,])
distances.eucl2 = rep(0,N)
for ( i in 1:N ) distances.eucl2[i] = eucl_dist(S2,u2[i,])
distances.eucl = c(distances.eucl1,distances.eucl2)
hist(distances.eucl, prob=TRUE)

sigma_eucl = sd(distances.eucl) / sqrt(e)
l_eucl = - log(2*pi) * e * N - log(sigma_eucl) * e * 2 * N - .5 * 2 * N * e

# X-axis grid
x = distances.eucl
x2 <- seq(min(x), max(x), length = 100)

# fun <- 2 * x2 * dchisq(x2^2, df=dim(u)[2])
# hist(distances.eucl, prob=TRUE, ylim= c(0,max(fun) ) )
# lines(x2, fun, col = 2, lwd = 2)

library("distory")
S1 = load_trees(paste(dir,"Species1.tre",sep=""))
u1 = load_trees(paste(dir,"Genes1_Depth0.7.nex",sep=""))
S2 = load_trees(paste(dir,"Species2.tre",sep=""))
u2 = load_trees(paste(dir,"Genes2_Depth0.7.nex",sep=""))

bhv_dist <- function(s,u){
  x = c(s,u)
  dist.multiPhylo(x)[1]
}

N = length(u1)
distances.bhv1 = rep(0,N)
for ( i in 1:N ) distances.bhv1[i] = bhv_dist(S1,u1[i])
distances.bhv2 = rep(0,N)
for ( i in 1:N ) distances.bhv2[i] = bhv_dist(S2,u2[i])
distances.bhv = c(distances.bhv1,distances.bhv2)
hist(distances.bhv)
hist(distances.bhv1)
hist(distances.bhv2)
