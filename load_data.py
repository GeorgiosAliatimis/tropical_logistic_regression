from Bio import Phylo

def load_dat(data_file_name):
    for u in Phylo.read(data_file_name,"nexus"):
        