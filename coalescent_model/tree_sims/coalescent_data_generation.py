from tree_sims import *
import random
import dendropy
import os 

rng = random.Random(1)
Rs = [.1,.2,.3,.4,.6,.8,1,1.2,1.4,1.6,1.8,2,4,6,8,10]
num_instances = 50
num_species_trees = 2 
num_gene_trees = 300 

dir_name = "coalescent_data"
if dir_name not in os.listdir("."):  
    os.mkdir(dir_name)

for R in Rs:
	print(R)
	if f"Depth{R}" not in os.listdir(dir_name):
		os.mkdir(f"{dir_name}/Depth{R}")

	for ind1 in range(num_instances):
		print(ind1)
		for ind2 in range(num_species_trees):
			sp_tree = species_tree(ntax=10,tree_depth = R,rng=rng)
			gene_trees = [gene_tree_from_species_tree(sp_tree) for _ in range(num_gene_trees)] 
			gene_trees = dendropy.TreeList(gene_trees)
			sp_tree = scale_tree(sp_tree)
			sp_tree.write(path=f"{dir_name}/Depth{R}/species_tree_{ind1+1}_{ind2+1}.nex",schema = "nexus")
			gene_trees.write(path=f"{dir_name}/Depth{R}/gene_trees_{ind1+1}_{ind2+1}.nex", schema = "nexus")
