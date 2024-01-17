rm(list=ls())
source("train_and_test.R")
source("load_data.R")
mode="classical" # tropical or classical for logistic regression
library(phangorn)

Rs = c(seq(0.2,2,0.2), seq(4,10,2))
e = 45
m = 10

Rs = c(0.1,1,10,0.2,0.4,0.6,0.8,1.2,1.4,1.6,1.8,2,4,6,8)

ROCs=list()
AUCs=list()
probs=list()
for (R in Rs){
  dir = paste("coalescent_model/tree_sims/data/Depth",R,sep="")
  reps = 1
  #reps = length(list.files(dir,patter="gene_trees")) %/% 2
  for(r in 1:reps){
    print(r)
    S1 = load_tree(paste(dir,"/species_tree_",2*r-1,".nex",sep=""))[[1]]
    S2 = load_tree(paste(dir,"/species_tree_",2*r,".nex",sep=""))[[1]]
    
    gene_trees_file1 = paste(dir,"/",list.files(dir,pattern=paste("gene_trees_",2*r-1,sep=""))[1],sep="")
    gene_trees_file2 = paste(dir,"/",list.files(dir,pattern=paste("gene_trees_",2*r,sep=""))[1],sep="")
    
    D1 = load_data(gene_trees_file1)
    D2 = load_data(gene_trees_file2)
    
    res <- train_and_test(D1,D2,mode=mode)
    ROCs[[paste(R)]] = res$ROC
    AUCs[[paste(R)]] = res$AUC
    probs[[paste(R)]] = res$probs
    omega_1 = res$pars[1:e]
    omega_2 = res$pars[(e+1):(2*e)]
    
    print(paste(R,res$AUC))
  }
}
