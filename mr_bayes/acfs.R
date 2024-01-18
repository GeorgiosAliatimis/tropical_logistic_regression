Ds = list()

samplefreq=16
ngen= 16e3
diagnfreq=4e3

# samplefreq=8
# ngen= 14e3
# diagnfreq=2e3

# samplefreq=4
# ngen= 15e3
# diagnfreq=1e3

N = round(ngen/diagnfreq)
for(i in 1:N) Ds[[i]] = read.csv(paste("primates_",ngen,"_",samplefreq,"_",diagnfreq,"/1/trees1/",i,".csv",sep="") ,sep=" ")
for(i in 1:N) names(Ds[[i]]) = 1:ncol(Ds[[1]])
D= Ds[[1]]
for(i in 2:N) D = rbind(D,Ds[[i]])

#D = D[500:nrow(D),]

LAG.MAX = 100
M = 1600

FUN = function (u) c(acf(u, plot = FALSE, lag.max=LAG.MAX)$acf)

a = matrix(0,nrow=ncol(D),ncol= LAG.MAX + 1)
for(i in 1:ncol(D)) a[i,] = FUN(D[,i]) 

plot((0:LAG.MAX) * samplefreq, colMeans(a), xlim=c(0,M ), ylim = c(0,1),ann=F,type='n')



for(k in 1:25){
    N = round(ngen/diagnfreq)
    for(i in 1:N) Ds[[i]] = read.csv(paste("primates_",ngen,"_",samplefreq,"_",diagnfreq,"/",k,"/trees1/",i,".csv",sep="") ,sep=" ")
    for(i in 1:N) names(Ds[[i]]) = 1:ncol(Ds[[1]])
    D= Ds[[1]]
    for(i in 2:N) D = rbind(D,Ds[[i]])

    FUN = function (u) c(acf(u, plot = FALSE, lag.max=LAG.MAX)$acf)

    a = matrix(0,nrow=ncol(D),ncol= LAG.MAX + 1)
    for(i in 1:ncol(D)) a[i,] = FUN(D[,i]) 

    lines((0:LAG.MAX) * samplefreq, colMeans(a), xlim=c(0,M ) )
}




samplefreq=8
ngen= 14e3
diagnfreq=2e3

LAG.MAX = 2* LAG.MAX
for(k in 1:25){
    N = round(ngen/diagnfreq)
    for(i in 1:N) Ds[[i]] = read.csv(paste("primates_",ngen,"_",samplefreq,"_",diagnfreq,"/",k,"/trees1/",i,".csv",sep="") ,sep=" ")
    for(i in 1:N) names(Ds[[i]]) = 1:ncol(Ds[[1]])
    D= Ds[[1]]
    for(i in 2:N) D = rbind(D,Ds[[i]])
    FUN = function (u) c(acf(u, plot = FALSE, lag.max=LAG.MAX)$acf)

    a = matrix(0,nrow=ncol(D),ncol= LAG.MAX + 1)
    for(i in 1:ncol(D)) a[i,] = FUN(D[,i]) 

    lines((0:LAG.MAX) * samplefreq, colMeans(a), xlim=c(0,M ) ,col="red")
}



samplefreq=4
ngen= 15e3
diagnfreq=1e3

lines((0:LAG.MAX) * samplefreq, colMeans(a), xlim=c(0,M ), ylim = c(0,1),ann=F,type='n')

LAG.MAX = 2* LAG.MAX
for(k in 1:25){
    N = round(ngen/diagnfreq)
    for(i in 1:N) Ds[[i]] = read.csv(paste("primates_",ngen,"_",samplefreq,"_",diagnfreq,"/",k,"/trees1/",i,".csv",sep="") ,sep=" ")
    for(i in 1:N) names(Ds[[i]]) = 1:ncol(Ds[[1]])
    D= Ds[[1]]
    for(i in 2:N) D = rbind(D,Ds[[i]])

    FUN = function (u) c(acf(u, plot = FALSE, lag.max=LAG.MAX)$acf)

    a = matrix(0,nrow=ncol(D),ncol= LAG.MAX + 1)
    for(i in 1:ncol(D)) a[i,] = FUN(D[,i]) 

    lines((0:LAG.MAX) * samplefreq, colMeans(a), xlim=c(0,M ), col="green")
}