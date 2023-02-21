###########################################
## Author :  Ruriko Yoshida
## Date   :  August 8th 2022
## Update :  August 9th 2022
## Program:  This code produces a tropical logistic regression for lungfish data
## Input  :  
## Output :  
## Execute: type in R as 
##
#############################################

# New library for calculating ROC curves
library(ROCR)
source("load_data.R")
source("test_model.R")
source("ulr.R")
source("generic_logistic_regression.R")

## NJ trees
## Reading lungfish data
n <- 10 ## number of leaves
T <- ape::read.tree("lungfish/lungfish_nj.txt") ## reading the input file of trees
N <- length(T) ## sample size of the data set
D <- matrix(rep(0, N*choose(n, 2)), N, choose(n, 2))
L <- T[[1]]$tip.label

for(i in 1:N){
    # u <- force.ultrametric(T[[i]], "nnls")
    u <- T[[i]]
    u$edge.length <- u$edge.length/max(u$edge.length)
    D0 <- cophenetic(u)[L, ]
    D[i, ] <- D0[lower.tri(t(D0))]
}

set.seed(876)




### NCUT with direct with L_2 metric

### Reading response variable
YY <- read.table("lungfish/fish_NJ_clustering_ncut_direct_cls2.txt", header = FALSE, sep = " ", dec = ".")

sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
D.train <- D[sampler,]
D.test <- D[-sampler,]

Y.train <- YY[sampler, 1]
Y.test <- YY[-sampler, 1]

res <- logistic(Y.train, D.train)
res

Y.hat <- Testing.Model(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC.NJ.direct <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
#plot(Logit.ROC.NJ.kpca, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC.NJ.direct <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
#Logit.AUC.NJ.kpca




### NCUT with TSNE with L_2 metric

### Reading response variable
YY <- read.table("lungfish/fish_NJ_clustering_ncut_tsne_D2_cls2.txt", header = FALSE, sep = " ", dec = ".")

sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
D.train <- D[sampler,]
D.test <- D[-sampler,]

Y.train <- YY[sampler, 1]
Y.test <- YY[-sampler, 1]

res <- logistic(Y.train, D.train)
res

Y.hat <- Testing.Model(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC.NJ.TSNE1 <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
#plot(Logit.ROC.NJ.kpca, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC.NJ.TSNE1 <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
#Logit.AUC.NJ.kpca




### NCUT with KPCA with L_2 metric

### Reading response variable
YY <- read.table("lungfish/fish_NJ_clustering_ncut_kpca_D2_cls2.txt", header = FALSE, sep = " ", dec = ".")

sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
D.train <- D[sampler,]
D.test <- D[-sampler,]

Y.train <- YY[sampler, 1]
Y.test <- YY[-sampler, 1]

res <- logistic(Y.train, D.train)
res

Y.hat <- Testing.Model(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC.NJ.kpca1 <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
#plot(Logit.ROC.NJ.kpca, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC.NJ.kpca1 <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
#Logit.AUC.NJ.kpca




### NCUT with TSNE with BHV metric

### Reading response variable
YY <- read.table("lungfish/fish_NJ_clustering_ncut_tsne_D3_cls2.txt", header = FALSE, sep = " ", dec = ".")

sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
D.train <- D[sampler,]
D.test <- D[-sampler,]

Y.train <- YY[sampler, 1]
Y.test <- YY[-sampler, 1]

res <- logistic(Y.train, D.train)
res

Y.hat <- Testing.Model(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC.NJ.TSNE2 <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
#plot(Logit.ROC.NJ.kpca, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC.NJ.TSNE2 <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
#Logit.AUC.NJ.kpca




### NCUT with KPCA with BHV metric

### Reading response variable
YY <- read.table("lungfish/fish_NJ_clustering_ncut_kpca_D3_cls2.txt", header = FALSE, sep = " ", dec = ".")

sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
D.train <- D[sampler,]
D.test <- D[-sampler,]

Y.train <- YY[sampler, 1]
Y.test <- YY[-sampler, 1]

res <- logistic(Y.train, D.train)
res

Y.hat <- Testing.Model(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC.NJ.kpca2 <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
#plot(Logit.ROC.NJ.kpca, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC.NJ.kpca2 <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
#Logit.AUC.NJ.kpca


#####  MLE tree

## Reading lungfish data
n <- 10 ## number of leaves
T <- ape::read.tree("lungfish/lungfish_mle.txt") ## reading the input file of trees
N <- length(T) ## sample size of the data set
D <- matrix(rep(0, N*choose(n, 2)), N, choose(n, 2))
L <- T[[1]]$tip.label

for(i in 1:N){
    # u <- force.ultrametric(T[[i]], "nnls")
    u <- T[[i]]
    u$edge.length <- u$edge.length/max(u$edge.length)
    D0 <- cophenetic(u)[L, ]
    D[i, ] <- D0[lower.tri(t(D0))]
}

set.seed(876)




### NCUT with direct with L_2 metric

### Reading response variable
YY <- read.table("lungfish/fish_MLE_clustering_ncut_direct_cls2.txt", header = FALSE, sep = " ", dec = ".")

sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
D.train <- D[sampler,]
D.test <- D[-sampler,]

Y.train <- YY[sampler, 1]
Y.test <- YY[-sampler, 1]

res <- logistic(Y.train, D.train)
res

Y.hat <- Testing.Model(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC.MLE.direct <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
#plot(Logit.ROC.MLE.kpca, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC.MLE.direct <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
#Logit.AUC.MLE.kpca




### NCUT with TSNE with L_2 metric

### Reading response variable
YY <- read.table("lungfish/fish_MLE_clustering_ncut_tsne_D2_cls2.txt", header = FALSE, sep = " ", dec = ".")

sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
D.train <- D[sampler,]
D.test <- D[-sampler,]

Y.train <- YY[sampler, 1]
Y.test <- YY[-sampler, 1]

res <- logistic(Y.train, D.train)
res

Y.hat <- Testing.Model(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC.MLE.TSNE1 <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
#plot(Logit.ROC.MLE.kpca, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC.MLE.TSNE1 <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
#Logit.AUC.MLE.kpca




### NCUT with KPCA with L_2 metric

### Reading response variable
YY <- read.table("lungfish/fish_MLE_clustering_ncut_kpca_D2_cls2.txt", header = FALSE, sep = " ", dec = ".")

sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
D.train <- D[sampler,]
D.test <- D[-sampler,]

Y.train <- YY[sampler, 1]
Y.test <- YY[-sampler, 1]

res <- logistic(Y.train, D.train)
res

Y.hat <- Testing.Model(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC.MLE.kpca1 <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
#plot(Logit.ROC.MLE.kpca, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC.MLE.kpca1 <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
#Logit.AUC.MLE.kpca




### NCUT with TSNE with BHV metric

### Reading response variable
YY <- read.table("lungfish/fish_MLE_clustering_ncut_tsne_D3_cls2.txt", header = FALSE, sep = " ", dec = ".")

sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
D.train <- D[sampler,]
D.test <- D[-sampler,]

Y.train <- YY[sampler, 1]
Y.test <- YY[-sampler, 1]

res <- logistic(Y.train, D.train)
res

Y.hat <- Testing.Model(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC.MLE.TSNE2 <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
#plot(Logit.ROC.MLE.kpca, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC.MLE.TSNE2 <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
#Logit.AUC.MLE.kpca




### NCUT with KPCA with BHV metric

### Reading response variable
YY <- read.table("lungfish/fish_MLE_clustering_ncut_kpca_D3_cls2.txt", header = FALSE, sep = " ", dec = ".")

sampler <- sample(nrow(D),trunc(nrow(D)*.80)) # samples index 
D.train <- D[sampler,]
D.test <- D[-sampler,]

Y.train <- YY[sampler, 1]
Y.test <- YY[-sampler, 1]

res <- logistic(Y.train, D.train)
res

Y.hat <- Testing.Model(res[[2]], D.test)
Y.class <- ifelse(Y.hat > 0.5, 1, 0)

# Calculate ROC statistics for our best logit model
Logit.ROC.MLE.kpca2 <- performance(prediction(Y.hat, Y.test), measure="tpr", x.measure="fpr")
#plot(Logit.ROC.MLE.kpca, lwd = 2, main = "ROC Curve for Logistic Regression Model")

# Calculate AUC - use @y.values to get AUC out
Logit.AUC.MLE.kpca2 <- performance(prediction(Y.hat, Y.test), 
                         measure="auc")@y.values
#Logit.AUC.MLE.kpca


### NJ
####  ROC
# Add direct
png("ROC_lungfish_NJ.png",width = 420, height = 420,)
plot(Logit.ROC.NJ.direct, lwd=2, main = "ROC for Tropical Logistic Regression with NJ")

# Add TNSE with L2
lines(attributes(Logit.ROC.NJ.TSNE1)$x.values[[1]], attributes(Logit.ROC.NJ.TSNE1)$y.values[[1]], 
      col="red", lwd=2)

# Add KPCA with L2
lines(attributes(Logit.ROC.NJ.kpca1)$x.values[[1]], attributes(Logit.ROC.NJ.kpca1)$y.values[[1]], 
      col="blue", lwd=2)

# Add TNSE with BHV
lines(attributes(Logit.ROC.NJ.TSNE2)$x.values[[1]], attributes(Logit.ROC.NJ.TSNE2)$y.values[[1]], 
      col="green", lwd=2)

# Add KPCA with BHV
lines(attributes(Logit.ROC.NJ.kpca2)$x.values[[1]], attributes(Logit.ROC.NJ.kpca2)$y.values[[1]], 
      col="brown", lwd=2)


#Add Legend
legend(x=.6, y=.4, c("Direct", "TNSE w/ L2", "KPCA w/ L2", "TNSE w/ BHV", "KPCA w/ BHV"), 
       col=c("black", "red", "blue", "green", "brown"), lwd=c(2,2,2,2,2,2))
dev.off()
#### AUC
Logit.AUC.NJ.direct
Logit.AUC.NJ.TSNE1
Logit.AUC.NJ.kpca1
Logit.AUC.NJ.TSNE2
Logit.AUC.NJ.kpca2

#### MLE
####  ROC
# Add direct
png("ROC_lungfish_MLE.png",width = 420, height = 420,)
plot(Logit.ROC.MLE.direct, lwd=2, main = "ROC for Tropical Logistic Regression with MLE")

# Add TNSE with L2
lines(attributes(Logit.ROC.MLE.TSNE1)$x.values[[1]], attributes(Logit.ROC.MLE.TSNE1)$y.values[[1]], 
      col="red", lwd=2)

# Add KPCA with L2
lines(attributes(Logit.ROC.MLE.kpca1)$x.values[[1]], attributes(Logit.ROC.MLE.kpca1)$y.values[[1]], 
      col="blue", lwd=2)

# Add TNSE with BHV
lines(attributes(Logit.ROC.MLE.TSNE2)$x.values[[1]], attributes(Logit.ROC.MLE.TSNE2)$y.values[[1]], 
      col="green", lwd=2)

# Add KPCA with BHV
lines(attributes(Logit.ROC.MLE.kpca2)$x.values[[1]], attributes(Logit.ROC.MLE.kpca2)$y.values[[1]], 
      col="brown", lwd=2)


#Add Legend
legend(x=.6, y=.4, c("Direct", "TNSE w/ L2", "KPCA w/ L2", "TNSE w/ BHV", "KPCA w/ BHV"), 
       col=c("black", "red", "blue", "green", "brown"), lwd=c(2,2,2,2,2,2))
dev.off()

#### AUC
Logit.AUC.MLE.direct
Logit.AUC.MLE.TSNE1
Logit.AUC.MLE.kpca1
Logit.AUC.MLE.TSNE2
Logit.AUC.MLE.kpca2

