"""
Auxiliary functions for the snakemake pipeline
"""
import pandas as pd

def get_introg_info(file,intro):
    reader = pd.read_csv(file,sep="\t",index_col=0)
    chr = reader.loc[intro]["Chromosome"]
    start = reader.loc[intro]["Start"]
    end = reader.loc[intro]["End"]
    return chr,start,end
