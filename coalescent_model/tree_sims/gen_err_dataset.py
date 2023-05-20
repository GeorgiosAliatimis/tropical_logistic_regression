from tree_sims import *
import random
import dendropy
import os 

#Using the same two species trees, we generate 1000 gene trees for each R. 


Rs = [.1,.2,.4,.6,.8,1,1.2,1.4,1.6,1.8,2,4,6,8,10]
num_gene_trees = 10000

dir_name = "gen_error_data"
if dir_name not in os.listdir("."):  
    os.mkdir(dir_name)

for R in Rs:
    print(R)
    rng = random.Random(1)
    sp_tree_1 = species_tree(ntax=10,tree_depth= R,rng=rng)
    sp_tree_2 = species_tree(ntax=10,tree_depth= R,rng=rng)
    
    for ind,sp_tree in enumerate([sp_tree_1,sp_tree_2]):
        gene_trees = [gene_tree_from_species_tree(sp_tree) for _ in range(num_gene_trees)] 
        gene_trees = dendropy.TreeList(gene_trees)
        gene_trees.write(path=f"{dir_name}/gene_trees_{R}_{ind+1}.nex", schema = "nexus")

scale_tree(sp_tree_1).write(path=f"{dir_name}/species_tree_1.nex",schema = "nexus")
scale_tree(sp_tree_2).write(path=f"{dir_name}/species_tree_2.nex",schema = "nexus")