rm(list=ls())
source("load_data.R")
source("train_and_test.R")

test_data <- function(file1,file2,model,color="black",is.plot=TRUE){
  res = train_and_test(file1,file2,model)
  # plot(res$ROC, lwd = 2, col=color,main = dir)
  if (is.plot){
    plot(res$ROC, lwd = 2, col=color,main = paste("ROC curves for R=",R,sep=""))
  } else {
    lines(res$ROC@x.values[[1]], res$ROC@y.values[[1]], lwd = 2, col=color)
  }
  res$AUC
}

Rs = c(0.1,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,4,6,8,10)
# Rs = c("025","05","1","2","5","10")
models = c("ulr","urat_lr","clr")
colors = c("green","red","black")
aucs = list()
for (R in Rs){
  dir = paste("./data/Depth",R,sep="")
  file1 = paste(dir,"/",list.files(dir,pattern="Genes1")[1],sep="")
  file2 = paste(dir,"/",list.files(dir,pattern="Genes2")[1],sep="")
  # file1 = paste("./data/coalescent_data/R",R,"genetrees_S1.dat",sep="")
  # file2 = paste("./data/coalescent_data/R",R,"genetrees_S2.dat",sep="")
  # png(paste("ROC_R=",R,".png"),         # File name
  #     width = 420, height = 420, # Width and height in inches
  #     bg = "white")
  for (model_ind in 1:3){
    model = models[model_ind] 
    color = colors[model_ind]
    auc = test_data(file1,file2,model,color,model_ind==1)
    aucs[[model]] = c(aucs[[model]], auc)
    print(paste("AUC=",auc))
  }
  # dev.off()
  print(paste("Validation of ",R,"complete."))
}

# save(aucs,file="aucs.Rdata")

inds= (Rs <= 2)
plot(Rs[inds],aucs[["clr"]][inds],ylab="AUC",xlab="R")
points(Rs[inds],aucs[["ulr"]][inds],col="green")
points(Rs[inds],aucs[["urat_lr"]][inds],col="red")

