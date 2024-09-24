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
common_file = out + "/" + prefix + "_common_introgressions_in_strains.csv"
common_file_nosub = out + "/" + prefix + "_common_introgressions_in_strains_nosub.csv"
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
    sp = os.path.basename(file).split("_")[0]
    intros = pd.read_csv(file, names=['SAPAorthID','SAPAcov','SAPAnorm','SACEorthID','SACEcov','SACEnorm','GO','Note',
              						  '%Identity','Gap_ID','SAPAstart','SAPAend','SAPAchr','SACEstart','SACEend','SACEchr',
              						  'Ratio','Heterozygosity','Subtelomeric','Has_ortholog'],skiprows=1,index_col=False)
    if del_subtelomeric:
        return set(intros.loc[intros['Subtelomeric'] == False]["SAPAorthID"])
    else:
        return set(intros["SAPAorthID"])

def sharedGenes(Dct,sp1,sp2,perc = False):
# This function takes a dictionary where the keys are the strains and the values are the sets of introgressed genes as obtained
# by the get_genes() function. {strain1:{gene1,gene2,...},strain2:{gene3,gene4,...},...}, it also requires the names of two strains.
# The function computes the intersection between the two sets of introgressed genes and outputs either the ammount of introgressed
# genes or the percentage. Depending of the value of the perc argument.
    genes1 = Dct[sp1]
    genes2 = Dct[sp2]
    noShared = len(genes1 & genes2)
    if perc:
        return noShared/len(genes1)
    else:
        return noShared

def subtelomeric_percentage(st,Dct_ns,Dct):
# This function takes a strain and the dictionaries with both all the genes and the nonsubtelomeric genes and computes the % of
# subtelomeric genes
    nosub = len(Dct_ns[st])
    total = len(Dct[st])
    subperc = (1-(nosub/total))*100
    return round(subperc,2)
# Open file with samples and set the Id's as index (this facilitates some computations)
sample_pd = pd.read_csv(samples).set_index("ID").fillna("NA")
# Obtain the strain names
strains = sample_pd.index
# Get all the introgressed genes files needed
genesFiles = get_files(folder,strains)

#Create the dictionary with the strains as keys and introgressed genes sets as values.
genesDict = {}
i = 0
for strain in strains:
    genesDict[strain] = get_genes(genesFiles[strain])
    # The if statement is to create the universe of introgresed genes
    if i == 0 :
        universe = genesDict[strain]
        i = 1
    else:
        universe = universe.union(genesDict[strain])
#Create the dictionary with the strains as keys and introgressed genes sets as values skipping subtelomeric genes
genesDict_nosub = {}
i = 0
for strain in strains:
    genesDict_nosub[strain] = get_genes(genesFiles[strain],del_subtelomeric=True)
    # The if statement is to create the universe of nonsubtelomeric introgresed genes
    if i == 0 :
        universe_nosub = genesDict_nosub[strain]
        i = 1
    else:
        universe_nosub = universe_nosub.union(genesDict_nosub[strain])

# Create data frames in which the data will go
head = sample_pd.columns.tolist()
common = pd.DataFrame(columns=head+list(universe),index=strains)
head_nosub = sample_pd.columns.tolist()
common_nosub = pd.DataFrame(columns=head_nosub+list(universe_nosub),index=strains)

for strain in strains:
    # Check if a certain gene is introgressed in all the strains
    for gene in universe:
        common["PhyloGroup_SACE469"][strain] = sample_pd["PhyloGroup_SACE469"][strain]
        present = int(gene in genesDict[strain])
        common[gene][strain] = present
    # Check if a certain nonsubtelomeric gene is introgressed in all the strains
    for gene in universe_nosub:
        common_nosub["PhyloGroup_SACE469"][strain] = sample_pd["PhyloGroup_SACE469"][strain]
        present = int(gene in genesDict_nosub[strain])
        common_nosub[gene][strain] = present

common.to_csv(common_file,header = True, sep = ',', index = True)
common_nosub.to_csv(common_file_nosub,header = True, sep = ',', index = True)
