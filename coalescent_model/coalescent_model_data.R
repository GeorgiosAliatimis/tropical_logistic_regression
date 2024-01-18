rm(list=ls())
source("train_and_test.R")
source("load_data.R")
mode="tropical" # tropical or classical for logistic regression
library(phangorn)

dir_name = "coalescent_data"
e = 45
m = 10

# number of times to run logistic regression on the same dataset, with different train/test splits
reps_per_example = 10 
Rs = sort(as.numeric(gsub("Depth","",list.files(dir_name))))
num_examples = length(list.files(paste(dir_name,"/Depth",min(Rs),sep=""),patter="gene_trees")) %/% 2

total_reps = reps_per_example * num_examples

AUCs=matrix(rep(0,total_reps * length(Rs)),nrow=length(Rs),ncol=total_reps)
for (i in 1:length(Rs)){
  R = Rs[i]
  dir = paste(dir_name,"/Depth",R,sep="")
  aucs = c()
  for(r in 1:num_examples){
    print(paste(R,r))
    S1 = load_tree(paste(dir,"/species_tree_",r,"_1.nex",sep=""))[[1]]
    S2 = load_tree(paste(dir,"/species_tree_",r,"_2.nex",sep=""))[[1]]
    gene_trees_file1 = paste(dir,"/gene_trees_",r,"_1.nex", sep="")
    gene_trees_file2 = paste(dir,"/gene_trees_",r,"_2.nex", sep="")
    
    D1 = load_data(gene_trees_file1)
    D2 = load_data(gene_trees_file2)
    
    aucs = c(aucs,sapply(1:reps_per_example,function(x) train_and_test(D1,D2,mode=mode)$AUC))
  }
  print(aucs)
  AUCs[i,] = aucs
}

library(MASS)
write.matrix(AUCs,file="aucs.csv",sep=",")
