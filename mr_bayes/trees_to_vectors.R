source("load_data.R")
library(MASS)
args = commandArgs(trailingOnly=TRUE)
dir = args[1]
input_file_1 = file.path(dir,"t1.nex")
input_file_2 = file.path(dir,"t2.nex")
output_file_1 = file.path(dir,"v1.csv")
output_file_2 = file.path(dir,"v2.csv")

D1 = load_data(input_file_1)
write.matrix(D1,file=paste(output_file_1,sep=""))
D2 = load_data(input_file_2)
write.matrix(D2,file=paste(output_file_2,sep=""))