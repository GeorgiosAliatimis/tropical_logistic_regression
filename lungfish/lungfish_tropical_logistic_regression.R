###########################################
## Author :  Ruriko Yoshida and Georgios Aliatimis
## Date   :  August 8th 2022
## Update :  April 4th 2023
## Program:  This code produces a tropical logistic regression for lungfish data
## Input  :  
## Output :  
## Execute: type in R as 
##
#############################################

# New library for calculating ROC curves
rm(list=ls())
library(ROCR)
source("load_data.R")
source("test_model.R")
source("lungfish/train_and_test.R")

## NJ trees
## Reading lungfish data

set.seed(12345)
dir = "lungfish/data/"

tree_files = c("lungfish_nj.txt","lungfish_mle.txt")
methods = c("NJ","MLE")

class_files = c("clustering_ncut_direct_cls2.txt",
                "clustering_ncut_tsne_D2_cls2.txt",
                "clustering_ncut_kpca_D2_cls2.txt",
                "clustering_ncut_tsne_D3_cls2.txt",
                "clustering_ncut_kpca_D3_cls2.txt")
classes = c("direct","TSNE1","kpca1","TSNE2","kpca2")

rocs = list()
aucs = list()
pars = list()
for(i in 1:length(tree_files)){
  tree_file = paste(dir,tree_files[i],sep="")
  method = methods[i]
  class_file_pref = paste(dir,"fish_",method,"_",sep="")
  for(j in 1:length(class_files)){
    class_file = paste(class_file_pref,class_files[i],sep="")
    class = classes[j]
    res = train_and_test(tree_file,class_file,model="fw")
    key = paste(method,".",class,sep="")
    rocs[[key]] = res$ROC
    aucs[key] = res$AUC
    pars[[key]] = res$pars
  }
}

### NJ
####  ROC
# Add direct
# png("ROC_lungfish_NJ.png",width = 420, height = 420,)
plot(rocs$NJ.direct, lwd=2, main = "ROC for Tropical Logistic Regression with NJ")

# Add TNSE with L2
lines(attributes(rocs$NJ.TSNE1)$x.values[[1]], attributes(rocs$NJ.TSNE1)$y.values[[1]], 
      col="red", lwd=2)

# Add KPCA with L2
lines(attributes(rocs$NJ.kpca1)$x.values[[1]], attributes(rocs$NJ.kpca1)$y.values[[1]], 
      col="blue", lwd=2)

# Add TNSE with BHV
lines(attributes(rocs$NJ.TSNE2)$x.values[[1]], attributes(rocs$NJ.TSNE2)$y.values[[1]], 
      col="green", lwd=2)

# Add KPCA with BHV
lines(attributes(rocs$NJ.kpca2)$x.values[[1]], attributes(rocs$NJ.kpca2)$y.values[[1]], 
      col="brown", lwd=2)


#Add Legend
legend(x=.6, y=.4, c("Direct", "TNSE w/ L2", "KPCA w/ L2", "TNSE w/ BHV", "KPCA w/ BHV"), 
       col=c("black", "red", "blue", "green", "brown"), lwd=c(2,2,2,2,2,2))
# dev.off()
#### AUC
print(aucs$NJ.direct)
print(aucs$NJ.TSNE1)
print(aucs$NJ.kpca1)
print(aucs$NJ.TSNE2)
print(aucs$NJ.kpca2)

#### MLE
####  ROC
# Add direct
# png("ROC_lungfish_MLE.png",width = 420, height = 420,)
plot(rocs$MLE.direct, lwd=2, main = "ROC for Tropical Logistic Regression with MLE")

# Add TNSE with L2
lines(attributes(rocs$MLE.TSNE1)$x.values[[1]], attributes(rocs$MLE.TSNE1)$y.values[[1]], 
      col="red", lwd=2)

# Add KPCA with L2
lines(attributes(rocs$MLE.kpca1)$x.values[[1]], attributes(rocs$MLE.kpca1)$y.values[[1]], 
      col="blue", lwd=2)

# Add TNSE with BHV
lines(attributes(rocs$MLE.TSNE2)$x.values[[1]], attributes(rocs$MLE.TSNE2)$y.values[[1]], 
      col="green", lwd=2)

# Add KPCA with BHV
lines(attributes(rocs$MLE.kpca2)$x.values[[1]], attributes(rocs$MLE.kpca2)$y.values[[1]], 
      col="brown", lwd=2)


#Add Legend
legend(x=.6, y=.4, c("Direct", "TNSE w/ L2", "KPCA w/ L2", "TNSE w/ BHV", "KPCA w/ BHV"), 
       col=c("black", "red", "blue", "green", "brown"), lwd=c(2,2,2,2,2,2))
# dev.off()

#### AUC
print(aucs$MLE.direct)
print(aucs$MLE.TSNE1)
print(aucs$MLE.kpca1)
print(aucs$MLE.TSNE2)
print(aucs$MLE.kpca2)


## Inferred trees 



labels = c("callorhinc","danio","gallus","homo",       "latimeria",  "leucoraja",  "lungfish",   "scyliorhin","takifugu", "xenopus"  )
plot_tree <- function(omega,title){
  dist_mat = matrix(0,m,m)
  dist_mat[lower.tri(dist_mat)]= omega - min(omega) + .1
  dist_mat = dist_mat + t(dist_mat)
  row.names(dist_mat) = labels
  colnames(dist_mat) = labels
  #u = phangorn::upgma(dist_mat)
  u = hclust(as.dist(dist_mat),method="complete")
  u = as.phylo(u)
  plot(u,main=title)
  u
}
m=10
e=45
titles = names(pars)
# tree_from_paper = ape::read.tree(text="(((((homo,gallus),xenopus), lungfish),latimeria), (danio, takifugu), ((scyliorhin, leucoraja), callorhinc));")
tree_from_paper = ape::read.tree(text="((((homo,gallus),xenopus), (lungfish,latimeria)), (danio, takifugu), ((scyliorhin, leucoraja), callorhinc));")

library(phangorn)
tree_dists = c()
for(i in 1:10){
  title = titles[i] 
  par = pars[i][[title]]
  omega_1 = par[1:e]
  # png(paste(title,"_1.png",sep=""))
  u1= plot_tree(omega_1,title)
  # dev.off()
  tree_dists = c(tree_dists,RF.dist(u1,tree_from_paper))
  
  if(length(par) > e+2){
    omega_2 = par[(e+1):(2*e)]
    # png(paste(title,"_2.png",sep=""))
    u2= plot_tree(omega_2,title)
    # dev.off()
    tree_dists = c(tree_dists,RF.dist(u2,tree_from_paper))
    # print(RF.dist(u1,u2))
  }
}
