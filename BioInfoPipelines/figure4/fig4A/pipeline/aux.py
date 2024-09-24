"""
Auxiliary functions for the snakemake pipeline
"""
import pandas as pd

def get_introgressions(file):
    reader = pd.read_csv(file,sep="\t")
    list_intros = list(reader["ID"])
    return list_intros

def get_introg_info(file,intro):
    reader = pd.read_csv(file,sep="\t",index_col=0)
    chr = reader.loc[intro]["Chromosome"]
    start = reader.loc[intro]["Start"]
    end = reader.loc[intro]["End"]
    return chr,start,end

def get_strains_share(file,intro):
    reader = pd.read_csv(file,sep="\t",index_col=False)
    strains = list(reader[intro].dropna())
    return strains
