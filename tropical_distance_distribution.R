source("load_data.R")
dir = "./data/Depth0.7/"

S = load_data(paste(dir,"Species1.tre",sep=""))
u = load_data(paste(dir,"Genes1_Depth0.7.nex",sep=""))

trop_dist <- function(x,y){
  max(x-y) - min(x-y)
}

eucl_dist <- function(x,y){
  sqrt( c(x-y) %*% c(x-y) )
}

N = dim(u)[1]
distances.trop = rep(0,N)
for ( i in 1:N ) distances.trop[i] = trop_dist(S,u[i,])
hist(distances.trop, prob=TRUE)

distances.eucl = rep(0,N)
for ( i in 1:N ) distances.eucl[i] = eucl_dist(S,u[i,])


# X-axis grid
x = distances.eucl
x2 <- seq(min(x), max(x), length = 100)

fun <- 2 * x2 * dchisq(x2^2, df=dim(u)[2])
hist(distances.eucl, prob=TRUE, ylim= c(0,max(fun) ) )
lines(x2, fun, col = 2, lwd = 2)
