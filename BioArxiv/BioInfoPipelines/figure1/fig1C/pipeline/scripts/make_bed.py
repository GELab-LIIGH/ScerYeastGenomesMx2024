####################################LIBRARIES###################################
import pandas as pd
import argparse, os
################################################################################
#################################PARSE ARGUMENTS################################
ag = argparse.ArgumentParser()
ag.add_argument("-g", "--gff", required = True, help = "Gff annotation file from S. cerevisiae")
ag.add_argument("-i", "--input", required = True, help = "List of orthologs to elaborate the tree with")
ag.add_argument("-o", "--out", required = True, help = "Output bed file name")
args = vars(ag.parse_args())
   ##################################INPUTS##################################
gff = args["gff"]
orth_list = args["input"]
   ##########################################################################
   #################################OUTPUTS##################################
bed_file = args["out"]
   ##########################################################################
################################################################################
####################################FUNCTIONS###################################
def gene_info(x):
    g_id = list(filter(lambda x: 'ID' in x,  x.split(";")))[0].split("=")[1]
    return g_id
################################################################################
if not os.path.isfile(bed_file):
    #Read gff
    gencode = pd.read_table(gff, comment="#",sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])
    #Extract important data
    gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)

    gencode_genes["ID"]  = gencode_genes.attribute.apply(lambda x: gene_info(x))
    gencode_genes = gencode_genes.sort_values(['seqname'], ascending=True).drop_duplicates('ID', keep='first').reset_index().drop('index', axis=1)
    #Read ortho list
    ortho = pd.read_table(orth_list, comment="#",sep = "\t", names = ['ID'])
    #Extract orthologs ffrom gff
    bed_pd = gencode_genes.merge(ortho)[['seqname', 'start', 'end']]
    #Write bed
    bed_pd.to_csv(bed_file, index=False, header = False, sep="\t")
