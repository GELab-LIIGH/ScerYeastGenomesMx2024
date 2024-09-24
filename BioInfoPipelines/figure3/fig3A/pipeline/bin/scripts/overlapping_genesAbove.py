import pandas as pd
import numpy as np
import pysam
import os
import subprocess
import argparse
import csv
#Also requires bedtools and tabix to be loaded

##Parse arguments
ag = argparse.ArgumentParser()
ag.add_argument("-t", "--txt", required = True, help = "txt file with genes")
ag.add_argument("-f", "--folder", required = True, help = "folder where the output files go")
ag.add_argument("-w", "--windows", required = True, help = "windows file to use")
ag.add_argument("-b", "--bamfile", required = True, help = "bam for coverage of genes")
args = vars(ag.parse_args())
txt_name = args["txt"]
out_folder = args["folder"]
win_csv = args["windows"]
bam = args["bamfile"]
##

ref = txt_name.split("/")[-1].split("_")[0]
gz = txt_name+".sorted.formated.bed.gz"
sample = win_csv.split("/")[-1].split("_")[0]
out_csv = out_folder + "/" + sample +"_"+ ref + "_genes_50aboveQm.csv"


def gencode_all_known_genes(a, tb):
    genes = []
    try:
        for region in tb.fetch(a['chr'], int(a['start']), int(a['end'])):
            if region:
                r = region.split('\t')
                overlap_len = overlap(int(a['start']), int(a['end']), int(r[1]), int(r[2]))
                ret_val = '{}'.format(r[3]) ### Percentage of the input interval that overlap with the gene
                genes.append(ret_val)
        if len(genes)>0:
            return ";".join(genes)
        else:
            return "NA"
    except ValueError:
        return "NA"

def overlap(q_st, q_end, res_st, res_end):
    o  = min(q_end, res_end)-max(q_st, res_st)
    return o

def fetch_gene_coords(g):
    if g in gencode_genes.index:
        return gencode_genes.loc[g]['seqname'], gencode_genes.loc[g]['start'], gencode_genes.loc[g]['end'],gencode_genes.loc[g]['Notes']   #gencode_genes.loc[g][['seqname', 'start', 'end']]
    else:
        return "NA", "NA", "NA","NA"

def retrieve_median(sample):
    file = os.path.dirname(win_csv)+"/"+sample+"_perReference.csv"
    with open(file,"r") as f:
        r = list(csv.reader(f,delimiter=","))
    hh = r.pop(0)
    dr = pd.DataFrame(r,columns=hh).sort_values(by="Median")
    return dr.Median[0]

df = pd.read_csv(win_csv, names=['chr', 'start', 'end', 'coverage'],skiprows=1)
#We only need the ref we are analyzing
df = df[df.chr.str.contains(ref)]
#windows with mean cov above 10 remain
df = df[df.coverage>1]

tbx = pysam.TabixFile(gz)
df['genes'] = df.apply(lambda x: gencode_all_known_genes(x[['chr', 'start', 'end']], tbx), axis=1)
#Remove non overlapping intervals
df = df[df['genes'] != "NA"].reset_index(drop=True)
if len(df) > 0:
    new_rows = []
    for i,r in df.iterrows():
        g_list = r['genes'].split(";")
        for g in g_list:
            g = g.replace(" ","")
            new_rows.append(np.append(r[['chr', 'start', 'end', 'coverage', 'genes']].values, g))

    df_perGene = pd.DataFrame(new_rows, columns=['chr', 'start', 'end','coverage', 'genes', 'ID']).reset_index().drop('index', axis=1)

    #No longer need those columns
    df_perGene = df_perGene.drop(["genes"], axis=1)
    ##Finding gene coordinates
    gencode_genes = pd.read_csv(txt_name,sep="\t",header=0,names=['seqname', 'start', 'end', 'attribute','ID','Notes']).set_index('ID')
    if len(df_perGene ) > 0:
        df_perGene['g_chr'], df_perGene['g_start'], df_perGene['g_end'], df_perGene['g_description'] = zip(*df_perGene['ID'].apply(lambda x: fetch_gene_coords(x)))
        df_perGene = df_perGene.drop_duplicates("ID",keep="first").reset_index().drop('index', axis=1)
        coverage = []
        above = []
        for i,gene in df_perGene.iterrows():
            chr = gene["g_chr"]
            start = gene["g_start"]
            end =gene["g_end"]
            get_cv = f"samtools depth -a -r {chr}:{start}-{end} {bam} | awk '{{sum+=$3}} END {{ print sum/NR}}'"
            try:
                ps_g = subprocess.Popen(get_cv,shell = True,stdout=subprocess.PIPE)
                ps_g.wait()
                cv = float(str(ps_g.communicate()[0]).split("'")[1].split("\\")[0])
            except:
                cv = 0
            coverage.append(cv)
            median = retrieve_median(sample)
            qmedian = int(median)*0.25
            get_abv = f"samtools depth -a -r {chr}:{start}-{end} {bam} | awk '$3 >= {qmedian} {{sum+=1}} END {{print sum/NR*100}}'"
            try:
                ps_a = subprocess.Popen(get_abv,shell = True,stdout=subprocess.PIPE)
                ps_a.wait()
                abv = float(str(ps_a.communicate()[0]).split("'")[1].split("\\")[0])
            except:
                abv = 0
            above.append(abv)
        df_perGene["coverage"] = coverage
        df_perGene["percAboveQMedian"] = above
        df_perGene = df_perGene[["ID","g_chr","g_start","g_end","coverage","percAboveQMedian","g_description"]]
        df_perGene = df_perGene[df_perGene.percAboveQMedian>50]
        if len(df_perGene) > 1 :
            df_perGene.to_csv(out_csv, index=False, header = True, sep=",")
