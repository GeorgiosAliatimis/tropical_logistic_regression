source("train_and_test.R")

args = commandArgs(trailingOnly=TRUE)
dir = args[1]
iter_file = args[2]
proportion_of_data = as.numeric(args[3])
output_file = file.path(dir,paste("aucs_",proportion_of_data,sep=""))

iterations = scan(iter_file)
aucs = c()
# p = as.numeric(args[2])
# N = as.numeric(args[3])
N = 10
thin <- function(D) D[sample(nrow(D),trunc(nrow(D)*p)),]

file1 = file.path(dir,"v1.csv") 
file2 = file.path(dir,"v2.csv")

D0 = as.matrix(read.csv(file1,sep=" ")) 
D1 = as.matrix(read.csv(file2,sep=" "))

aucs = scan(output_file)
for(end in iterations){
	print(end)
	if(end < 2e4) next
	# end = i * nrow(D0)%/%num_generations 
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
	start = end + 1
	print(paste(end,mean(auc),"+-", sd(auc)/sqrt(N)))
}
print(aucs)
write.table(aucs,file=output_file,col.names=F,row.names=F)