import pyreadr
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.metrics import silhouette_score
import sys

# Get matrix D
D = pyreadr.read_r("./mle_gene_trees")
D = D["D"]
D = np.array(D)

clust = DBSCAN(eps = .6, min_samples=8).fit(D)
labs = -clust.labels_
print(sum(labs == 0))
print(sum(labs == 1))
print(sum(labs < 0))
print(silhouette_score(D,clust.labels_))
labs = [int(lab != 0) for lab in labs]

lines = [f"{lab} {1-lab}\n" for lab in labs]

file = open("./lungfish/data/fish_MLE_dbscan_clustering.txt","w")
file.writelines(lines)
file.close()

import random
file = open("./lungfish/data/fish_MLE_random_clustering.txt","w")
labs = [random.randint(0,1) for _ in range(len(labs))]
lines = [f"{lab} {1-lab}\n" for lab in labs]
file.writelines(lines)
file.close()
