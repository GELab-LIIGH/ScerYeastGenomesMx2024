import pandas as pd
import os
import subprocess
import argparse
#Also requires bedtools and tabix to be loaded

##Parse arguments
ag = argparse.ArgumentParser()
ag.add_argument("-g", "--gff", required = True, help = "gff file of the reference genome")
ag.add_argument("-f", "--folder", required = True, help = "folder where the output file goes")
args = vars(ag.parse_args())
gff = args["gff"]
out_folder = args["folder"]
##

ref = gff.split("/")[-1].split("_")[0]
txt_name = out_folder +os.path.basename(gff).split(".gff")[0] + ".txt"
gz = txt_name+".sorted.formated.bed.gz"

def gene_info(x):
    try:
        g_id = list(filter(lambda x: 'ID' in x,  x.split(";")))[0].split("=")[1]
    except:
        g_id = list(filter(lambda x: 'gene_id' in x,  x.split(";")))[0].split("=")[1]
    try:
        g_note = list(filter(lambda x: 'Note' in x,  x.split(";")))[0].split("=")[1]
    except:
        try:
            g_note = list(filter(lambda x: 'gene_biotype' in x,  x.split(";")))[0].split("=")[1]
        except:
            try:
                g_note = list(filter(lambda x: 'biotype' in x,  x.split(";")))[0].split("=")[1]
            except:
                g_note = "NA"
    return (g_id,g_note)

gencode = pd.read_table(gff, comment="#",sep = " ", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])
if pd.isna(gencode.source[0]):
    gencode = pd.read_table(gff, comment="#",sep = "\t", names = ['seqname', 'source', 'feature', 'start' , 'end', 'score', 'strand', 'frame', 'attribute'])

gencode_genes = gencode[(gencode.feature == "gene")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)

if len(gencode_genes) == 0:
    gencode_genes = gencode[(gencode.feature == "CDS")][['seqname', 'start', 'end', 'attribute']].copy().reset_index().drop('index', axis=1)

gencode_genes["ID"],gencode_genes["Notes"]  = zip(*gencode_genes.attribute.apply(lambda x: gene_info(x)))
gencode_genes = gencode_genes.sort_values(['seqname'], ascending=True).drop_duplicates('ID', keep='first').reset_index().drop('index', axis=1)
if not os.path.isfile(txt_name):
    gencode_genes.to_csv(txt_name, index=False, header = False, sep="\t")
bedfile = txt_name+".sorted.formated.bed"

if not os.path.isfile(bedfile+".gz"):
    mk_bed = f"cut -f 1,2,3,5 {txt_name} | sortBed -i > {bedfile}"
    ps = subprocess.Popen(mk_bed,shell=True)
    ps.wait()
    mk_bgzp = ["bgzip",bedfile]
    subprocess.call(mk_bgzp)
    mk_tbx = ["tabix","-p","bed",gz]
    subprocess.call(mk_tbx)
