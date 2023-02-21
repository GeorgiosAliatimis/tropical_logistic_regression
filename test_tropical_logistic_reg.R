###########################################
## Author :  Ruriko Yoshida
## Date   :  April 23rd 2022
## Update :  July 23th 2022
## Program:  This code produces a tropical logistic regression
## Input  :  
## Output :  
## Execute: type in R as 
##
#############################################

#source("tropical_LM.R")
# source("tropical_logistic_regression.R")
library(ape)
library(phangorn)
library(phytools)
source("tropical_PCA.R")

# New library for calculating ROC curves
library(ROCR)

Testing.tropical.LM <- function(pars, D){
    ## omega is the normal vector for the best-fit Stiefel tropical hyperplane
    ## D is a nxe matrix.  the ith row is the vector for the predictors for the ith obs.
    d <- dim(D)
    Y.hat <- rep(0, d[1])
    D = D - rowMeans(D)   
    for(i in 1:d[1])
        Y.hat[i] <- prob(pars, D[i,])
    return(Y.hat)
}

train.size <- 800
test.size <- 200

set.seed(876)
### R = 1
# T1 <- ape::read.nexus("data/taxa6R1genetrees_S1.dat")
# T2 <- ape::read.nexus("data/taxa6R1genetrees_S2.dat")
R = "2"
T1 <- ape::read.nexus(paste("data/R",R,"genetrees_S1.dat",sep=""))
T2 <- ape::read.nexus(paste("data/R",R,"genetrees_S2.dat",sep=""))
n <- 10 ## number of leaves

L <- T2[[1]]$tip.label
d <- choose(n, 2)
N1 <- length(T1) ## sample size of the data set
D1 <- matrix(rep(0, N1*choose(n, 2)), N1, choose(n, 2))
N2 <- length(T2) ## sample size of the data set
D2 <- matrix(rep(0, N2*choose(n, 2)), N2, choose(n, 2))

for(i in 1:N1){
    u <- force.ultrametric(T1[[i]], "nnls")
    u$edge.length <- u$edge.length/max(u$edge.length)
    D0 <- cophenetic(u)[L, ]
    D1[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}
for(i in 1:N2){
    u <- force.ultrametric(T2[[i]], "nnls")
    u$edge.length <- u$edge.length/max(u$edge.length)
    D0 <- cophenetic(u)[L, ]
    D2[i, ] <- normaliz.tree(D0[lower.tri(t(D0))], 2)
}

## training set
DD1 <- D1[1:train.size, ]
DD2 <- D2[1:train.size, ]

## test set
DDD1 <- D1[(train.size+1):(train.size+test.size), ]
DDD2 <- D2[(train.size+1):(train.size+test.size), ]

D.train <- rbind(DD1, DD2)
Y.train <- c(rep(0, train.size), rep(1, train.size))

D.test <- rbind(DDD1, DDD2)
Y.test <- c(rep(0, test.size), rep(1, test.size))

source("tropical_logistic_regression_alpha.R")
set.seed(1)
res <- tropical.logistic(Y.train, D.train)
res

Y.hat <- Testing.tropical.LM(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
plot(Logit.ROC, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
Logit.AUC
