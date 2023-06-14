library(ROCR)
source("load_data.R")
source("test_model.R")

train_and_test <- function(tree_file,class_file,model="fw",method="two_species"){
  source(paste("methods/",model,".R",sep=""))
  
  D = load_data(tree_file)
  Y = read.table(class_file, header = FALSE, sep = " ", dec = ".")
  
  sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
  D.train <- D[sampler,]
  D.test <- D[-sampler,]
  
  Y.train <- Y[sampler, 1]
  Y.test <- Y[-sampler, 1]
  
  res <- logistic(Y.train, D.train,method=method)
  
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
