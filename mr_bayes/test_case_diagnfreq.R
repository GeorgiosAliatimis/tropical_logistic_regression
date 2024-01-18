source("train_and_test.R")

args = commandArgs(trailingOnly=TRUE)
dir = "primates_4_4000"

GENERATIONS=length(list.files(file.path(dir,"trees1")))

p = 1
N = 10
thin <- function(D) D[sample(nrow(D),trunc(nrow(D)*p)),]

i=3
treefile=paste(i,".csv",sep="")
file1 = file.path(dir,"trees1",treefile) 
file2 = file.path(dir,"trees2",treefile) 

D0 = as.matrix(read.csv(file1,sep=" "))
D1 = as.matrix(read.csv(file2,sep=" "))

i=4
treefile=paste(i,".csv",sep="")
file1 = file.path(dir,"trees1",treefile) 
file2 = file.path(dir,"trees2",treefile) 

names(D0) = 1:ncol(D0)
names(D1) = 1:ncol(D1)
D0 = rbind(D0, as.matrix(read.csv(file1,sep=" ")))
D1 = rbind(D1, as.matrix(read.csv(file2,sep=" ")))


aucsk = c()
for( k in 7:15){
	aucs = c()
	for(i in 1:3){
		D0i = D0[(1+ (i-1) * nrow(D0) %/% k ):(i* nrow(D0) %/% k),]
		D1i = D1[(1+ (i-1) * nrow(D1) %/% k ):(i* nrow(D1) %/% k),]
		auc = c()
		for ( j in 1:N ){
			out = train_and_test(thin(D0i),thin(D1i))
			auc = c(auc,out$AUC[[1]])
		}
		aucs = c(aucs, mean(auc))
	}
	print(aucs)
	aucsk = c(aucsk,mean(aucs))
}


print(i)
print(mean(auc))
print(sd(auc)/sqrt(N))

print(aucs)
output_file = file.path(dir,"aucs_tlr")
write.table(aucs,file=output_file,col.names=F,row.names=F)