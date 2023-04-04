###########################################
## Author :  Ruriko Yoshida
## Date   :  November 25th 2022
## Update :  November 25th 2022
## Program:  This code produces a tropical logistic regression
## Input  :  
## Output :  
## Execute: type in R as 
##
#############################################

# New library for calculating ROC curves
rm(list=ls())
library(ROCR)
source("load_data.R")
source("ulr.R")
source("test_model.R")

# Loading data  
R="025"
datafile0 = paste("./data/coalescent_data/R",R,"genetrees_S1.dat",sep="")
datafile1 = paste("./data/coalescent_data/R",R,"genetrees_S2.dat",sep="")

D0 = load_data(datafile0)
D1 = load_data(datafile1)

stopifnot(dim(D0) == dim(D1))

# 80-20 validation split
N = dim(D0)[1]
train.size <- floor(0.8 * N)
test.size <- N - train.size

## (D,Y) for training and testing
train_set = 1:train.size
test_set  = (train.size+1):(train.size+test.size) 

D.train <- rbind(D0[train_set,], D1[train_set,])
Y.train <- c(rep(0, train.size), rep(1, train.size))

D.test <- rbind(D0[test_set,], D1[test_set,])
Y.test <- c(rep(0, test.size), rep(1, test.size))

D.train <- normalize_data(D.train)
D.test  <- normalize_data(D.test) 
set.seed(1) 
res <- logistic(Y.train,D.train,logistic_gradient, pars_gen,I)

pars <- res$omega
logit.cost.function(pars,Y.train,D.train)

Y.hat <- Testing.Model(pars, D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)
# Calculate ROC statistics for our best logit model
Logit.ROC <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
plot(Logit.ROC, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
print(Logit.AUC)