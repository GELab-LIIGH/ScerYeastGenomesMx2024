import os
import argparse
import pandas as pd

#Parse arguments
ag = argparse.ArgumentParser()
ag.add_argument("-s", "--samples", required = True, help = """CSV file with the data of the samples to work.\n
The document must contain the location or phylogenetic groups information to group the strains""")
ag.add_argument("-o", "--out", required = True, help = "Folder where the output goes")
ag.add_argument("-f", "--folder", required = True, help = "Folder where the introgressed genes files is")
ag.add_argument("-p", "--prefix", required = True, help = """Prefix which wil be added to the output files, allowing multiple matrices.\n
This script generates 9 ouput files, the 'nosub.csv'files intentionally ignore subtelomeric genes.""")

args = vars(ag.parse_args())
samples = args["samples"]
out = args["out"]
folder = args["folder"]
prefix = args["prefix"]

############################ NAMES OF OUTPUT FILES ############################
absolutes_file = out + "/" + prefix + "_absolutes.csv"
bygroup_file = out + "/" + prefix + "_frequencies.csv"
###############################################################################

def get_files(path,samples):
# This function takes the path where the files with the introgressions are and the list of samples to
# retrieve the files with the introgressions
# This function has as output a dictionary where the keys are the strains and the values are the introgression
# file locarions: {strain1:'/location/strain1_introSAPA_All_noOrth.csv',strain2:'/location/strain2_introSAPA_All_noOrth.csv',...}
    dirs = os.listdir(path)
    files = {}
    for entry in dirs:
        full_path = os.path.join(path, entry)
        if os.path.isdir(full_path):
            pass
        else:
            if "_introSAPA_All_noOrth.csv" in full_path :
                files[entry.split("_")[0]] = full_path
    return files

def get_genes(file,del_subtelomeric = False):
# This function takes the introgression file and outputs the introgessed genes as a set.
# Output example: {gene1, gene2, gene3, ...}
# Set subtelomeric = True to erase all subtelomeric genes.
    #Read introgressions file
    intros = pd.read_csv(file, names=['SAPAorthID','SAPAcov','SAPAnorm','SACEorthID','SACEcov','SACEnorm','GO','Note',
              						  '%Identity','Gap_ID','SAPAstart','SAPAend','SAPAchr','SACEstart','SACEend','SACEchr',
              						  'Ratio','Heterozygosity','Subtelomeric','Has_ortholog'],skiprows=1,index_col=False)

    #We only need the name and the heterozigosity
    if del_subtelomeric:
        #Filter out subtelomeric
        intros_het = intros.loc[intros['Subtelomeric'] == False][["SAPAorthID",'Heterozygosity']]
    else:
        intros_het = intros[["SAPAorthID",'Heterozygosity']]
    #This dict lokks something like:
    #{introgresion1:heterozigosity1,introgresion2:heterozigosity2,...,introgresionN:heterozigosityN}
    intros_dict = {i["SAPAorthID"]:i["Heterozygosity"] for n, i in intros_het.iterrows()}
    return intros_dict

# Open file with samples and set the Id's as index (this facilitates some computations)
sample_pd = pd.read_csv(samples).set_index("ID").fillna("NA")
# Obtain the strain names
strains = sample_pd.index
# Get all the introgressed genes files needed
genesFiles = get_files(folder,strains)

#Create the dictionary with the strains as keys and as values, dictionaries with the gene as
#keys and its heterozigosity as value
genesDict = {}
#Create the dictionary were the universe of introgressed genes will go
universe = set()
for strain in strains:
    genesDict[strain] = get_genes(genesFiles[strain])
    universe = universe.union(set(genesDict[strain].keys()))

# Create data frames in which the data will go
groups1 = sample_pd.PhyloGroup_SACE469.unique().tolist()
bygroup = pd.DataFrame(columns=groups1, index=universe).fillna(0)

#Obtain counts of each gene in a specific population
for strain in strains:
    # Get groups of the strain
    g1 = sample_pd.loc[strain]["PhyloGroup_SACE469"]
    # Check if a certain gene is introgressed in all the strains
    for gene in universe:
        present = int(gene in genesDict[strain].keys())
        # If the gene is present add 1 to the ammount of strains of that group that share the introgression
        bygroup[g1][gene] += present

#Obtain size of populations
pop_s = sample_pd["PhyloGroup_SACE469"].value_counts()

rel = bygroup.copy()

#Obtain frequency of each gene in a specific population
for g in groups1:
    rel[g] = rel[g]/pop_s[g]

bygroup.to_csv(absolutes_file,header = True, sep = ',', index = True)
rel.to_csv(bygroup_file,header = True, sep = ',', index = True)
