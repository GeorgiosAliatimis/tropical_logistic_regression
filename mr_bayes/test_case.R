source("train_and_test.R")

args = commandArgs(trailingOnly=TRUE)
dir = args[1]
iter_file = args[2]
proportion_of_data = as.numeric(args[3])

input_file_1 = file.path(dir,"v1.csv")
input_file_2 = file.path(dir,"v2.csv")
output_file = file.path(dir,"aucs")

iterations = scan(iter_file)
aucs = c()
N = 10
thin <- function(D) D[sample(nrow(D),trunc(nrow(D)*p)),]

D0 = as.matrix(read.csv(input_file_1,sep=" ")) 
D1 = as.matrix(read.csv(input_file_2,sep=" "))

for(end in iterations){
	print(end)
	start = round( (1-proportion_of_data) * end)
	indices = start:end
	D0i = D0[indices,]
	D1i = D1[indices,]
	p = min(1, 300/nrow(D0i))
	auc = c()
	for ( j in 1:N ){
		out = train_and_test(thin(D0i),thin(D1i))
		auc = c(auc,out$AUC[[1]])
	}
	aucs = c(aucs, mean(auc))
	print(paste(end,mean(auc),"+-", sd(auc)/sqrt(N)))
}
print(aucs)
write.table(aucs,file=output_file,col.names=F,row.names=F)
