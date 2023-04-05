rm(list=ls())
source("train_and_test.R")
model="ulr"
Rs = c(0.1,0.6,0.7,0.8,0.9,1,1.2,1.4,1.6,1.8,2,4,6,8,10)
c1 <- rgb(173,216,230,max = 255, alpha = 80, names = "lt.blue")
c2 <- rgb(255,192,203, max = 255, alpha = 80, names = "lt.pink")
for(R in Rs){
  dir = paste("./data/Depth",R,sep="")
  file1 = paste(dir,"/",list.files(dir,pattern="Genes1")[1],sep="")
  file2 = paste(dir,"/",list.files(dir,pattern="Genes2")[1],sep="")
  res = train_and_test(file1,file2,model)
  P = res$probs
  N = length(P)/2
  
  print( ( sum(P[P<1/2]) + sum(1-P[P>1/2]) ) / length(P) )
  breaks = seq(0,1,.05)
  hgA <- hist(P[1:N],  plot = FALSE,breaks=breaks) # Save first histogram data
  hgB <- hist(P[(N+1):(2*N)], plot = FALSE,breaks=breaks) # Save 2nd histogram data
  
  png(file=paste("hist_estimated_prob_R=",R,".png"))
  plot(hgA, col = c1,main=paste("R=",R),
       xlab="Estimated probability p",
       xlim=c(min(P),max(P)),
       ylim= c(0,max(hgA$counts,hgB$counts))) # Plot 1st histogram using a transparent color
  plot(hgB, col = c2, add = TRUE) # Add 2nd histogram using different color
  dev.off()
}