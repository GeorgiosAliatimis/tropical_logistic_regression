rm(list=ls())
source("coalescent_model/train_and_test.R")
source("load_data.R")
library(phangorn)
model="ulr"
Rs = c(seq(0.2,2,0.2), seq(4,10,2))
e = 45
m = 10

plot_species_and_inferred_tree <- function(omega,st,labels,title_spe="",title_inf=""){
  # Species tree
  dist_mat <- cophenetic(st)
  dist_mat <- dist_mat[labels,labels]
  dist_mat <- dist_mat / max(dist_mat)
  s_vec <- dist_mat[lower.tri(dist_mat)]
  s = phangorn::upgma(dist_mat)
  # s = hclust(as.dist(dist_mat),method="complete")
  # s = as.phylo(s)
  plot(s,main=title_spe)
  ms = min(dist_mat[dist_mat>0])
  # Inferred tree
  dist_mat = matrix(0,m,m)
  dist_mat[lower.tri(dist_mat)]= omega - min(omega) + ms
  dist_mat = dist_mat + t(dist_mat)
  row.names(dist_mat) = labels
  colnames(dist_mat) = labels
  u = phangorn::upgma(dist_mat)
  # u = hclust(as.dist(dist_mat),method="complete")
  # u = as.phylo(u)
  plot(u,main=title_inf)
  # all.equal.phylo(s,u,use.edge.length=F)
  #summary(lm(omega ~ s_vec))$r.squared
  if (all.equal.phylo(s,u,use.edge.length=F)){
    print(paste("Topology of ",title_spe,"matches that of ", title_inf))
  } else {
    # print(paste("Topology of ",title_spe,"does not match that of ", title_inf))
    print(paste("Tree distance between ", title_spe," and ", title_inf, " is ", RF.dist(s,u)  ))
  }
  c(RF.dist(s,u),trop_dist(s_vec,omega))
}

tip_labels = c("a","b","c","d","e","f","g","h","i","j")

Rs = c(0.1,1,10,0.2,0.4,0.6,0.8,1.2,1.4,1.6,1.8,2,4,6,8)

tree_dists.RF = list()
tree_dists.trop = list()
ROCs=list()
AUCs=list()
probs=list()
for (R in Rs){
  dir = paste("coalescent_model/tree_sims/data/Depth",R,sep="")
  reps = 1
  #reps = length(list.files(dir,patter="gene_trees")) %/% 2
  tree_dists.RF[[paste(R)]] = rep(0,2*reps)
  tree_dists.trop[[paste(R)]] = rep(0,2*reps)
  for(r in 1:reps){
    print(r)
    S1 = load_tree(paste(dir,"/species_tree_",2*r-1,".nex",sep=""))[[1]]
    S2 = load_tree(paste(dir,"/species_tree_",2*r,".nex",sep=""))[[1]]
    
    gene_trees_file1 = paste(dir,"/",list.files(dir,pattern=paste("gene_trees_",2*r-1,sep=""))[1],sep="")
    gene_trees_file2 = paste(dir,"/",list.files(dir,pattern=paste("gene_trees_",2*r,sep=""))[1],sep="")
    
    res <- train_and_test(gene_trees_file1,gene_trees_file2, model = "fw")
    print(res$log_lik_val)
    ROCs[[paste(R)]] = res$ROC
    AUCs[[paste(R)]] = res$AUC
    probs[[paste(R)]] = res$probs
    omega_1 = res$pars[1:e]
    omega_2 = res$pars[(e+1):(2*e)]
    
    d1 <- plot_species_and_inferred_tree(omega_1,S1,tip_labels,paste("Species tree 1, R=",R), paste("Inferred tree 1, R=",R))
    d2 <- plot_species_and_inferred_tree(omega_2,S2,tip_labels,paste("Species tree 2, R=",R), paste("Inferred tree 2, R=",R))
    
    tree_dists.RF[[paste(R)]][2*r-1] = d1[1]
    tree_dists.RF[[paste(R)]][2*r] = d2[1]
    tree_dists.trop[[paste(R)]][2*r-1] = d1[2]
    tree_dists.trop[[paste(R)]][2*r] = d2[2]
  }
}
# save(list = c("Rs","ROCs","AUCs","tree_dists.RF","tree_dists.trop"),file="python_fw_all_Rs_big.RData")
