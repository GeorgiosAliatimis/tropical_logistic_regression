library(ROCR)
source("load_data.R")
source("test_model.R")

train_and_test <- function(datafile0,datafile1,model="ulr",method="two_species"){
  source(paste("methods/",model,".R",sep=""))
  
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
  
  set.seed(1) 
  res <- logistic(Y.train,D.train,method=method)
  
  pars <- res$omega

  Y.hat <- Testing.Model(pars, D.test)
  Y.class <- ifelse(Y.hat > 0.5, 1, 0)
  
  Logit.ROC <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
  # plot(Logit.ROC, lwd = 2, main = "ROC Curve for Logistic Regression Model")
  
  # Calculate AUC - use @y.values to get AUC out
  Logit.AUC <- performance(prediction(Y.hat, Y.test), 
                           measure="auc")@y.values
  list(pars=pars,probs=Y.hat,ROC=Logit.ROC,AUC=Logit.AUC)
}
