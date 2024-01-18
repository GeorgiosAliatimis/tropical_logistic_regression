source("train_and_test.R")

args = commandArgs(trailingOnly=TRUE)
dir = args[1]

GENERATIONS=length(list.files(file.path(dir,"trees1")))

aucs = c()

p= as.numeric(args[2])
N = as.numeric(args[3])
output_file = args[2]
thin <- function(D) D[sample(nrow(D),trunc(nrow(D)*p)),]
for(i in 1:4){
	treefile=paste(i,".csv",sep="")
	file1 = file.path(dir,"trees1",treefile) 
	file2 = file.path(dir,"trees2",treefile) 

	if(i==1) D0 = as.matrix(read.csv(file1,sep=" "))
    else D0 = rbind(D0,as.matrix(read.csv(file1,sep=" ")))
    if(i==1) D1 = as.matrix(read.csv(file2,sep=" "))
    else D1 = rbind(D1,as.matrix(read.csv(file2,sep=" ")))
}
print(dim(D0));print(dim(D1))
auc = c()
for ( j in 1:N ){
    out = train_and_test(thin(D0),thin(D1))
    auc = c(auc,out$AUC[[1]])
}
aucs = c(aucs, mean(auc))
print(i)
print(mean(auc))
print(sd(auc)/sqrt(N))

print(aucs)
