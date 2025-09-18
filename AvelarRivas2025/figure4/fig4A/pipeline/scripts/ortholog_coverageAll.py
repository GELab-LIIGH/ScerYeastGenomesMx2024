from __future__ import print_function, division
import pandas as pd
import numpy as np
import pysam, os, csv, argparse, subprocess, sys
try:
    from Bio import pairwise2
    from Bio import SeqIO
    from Bio.SubsMat import MatrixInfo as matlist
except ImportError as exception:
    print("[!] Could not import Biopython modules", file=sys.stderr)
    raise exception

##Parse arguments
ag = argparse.ArgumentParser()
ag.add_argument("-f", "--folder", required = True, help = "folder with files")
ag.add_argument("-a", "--sapa", required = True, help = "sapa genes with coverage on sample")
ag.add_argument("-e", "--sace", required = True, help = "sace genes with coverage in sample")
ag.add_argument("-b", "--bam", required = True, help = "bam file against the concatenated reference")
ag.add_argument("-o", "--orths", required = True, help = "csv with samples")
ag.add_argument("-O", "--out_folder", required = True, help = "folder where the output file will go")
ag.add_argument("-s", "--sample", required = True, help = "sample")
args = vars(ag.parse_args())

folder = args["folder"]
sapafile = args["sapa"]
sacefile = args["sace"]
bam = args["bam"]
orths = args["orths"]
out_folder = args["out_folder"]
sample = args["sample"]
#
sapatxt = folder+"SAPA_YPS138_v1.txt"
sacetxt = folder+"SACE_S288C_v1_allChr.txt"
fastaSACE = folder+"SACE_S288C_v1_allChr.fna"
fastaSAPA = folder+"SAPA_YPS138_v1.fna"
subtelomeric_regions = folder+"SAPA_YPS138_v1_subtelomeric_regions.csv"
#

## Functions
#Get gene ontology # Edited to obtain any attribute, not just GO.
def gene_info(x,attr):
    g_ont = list(filter(lambda x: attr in x,  x.split(";")))[0].split("=")[1]
    return g_ont
def subtelomeric(gchr,gstart,gend):
	lh_end = subt_r["L_end"][gchr]
	rh_start = subt_r["R_end"][gchr]
	if int(gstart) > lh_end and int(gend) < rh_start:
		return False
	else:
		return True

#Reads csv and returns 2 dimensional array
def csv2list(path):
    file = open(path,"r")
    i = 0
    arr = []
    for line in file:
        if  i == 0:
            pass
            i+=1
        else:
            l = line.strip().split(',')
            arr.append(l)
    file.close()
    return arr
#Find orthologs from sace and sapa and extract data from txt's
def find_orthologs(sapaf,sacef):
	seqIterSACEtmp = SeqIO.parse(fastaSACE, 'fasta')
	seqIterSAPAtmp = SeqIO.parse(fastaSAPA, 'fasta')
	seqIterSACE = SeqIO.to_dict(seqIterSACEtmp)
	seqIterSAPA = SeqIO.to_dict(seqIterSAPAtmp)
	prevKeys = seqIterSAPA.keys()
	for k in prevKeys:
		seqIterSAPA[k.split("|")[0]] = seqIterSAPA.pop(k)
	orthologs = pd.DataFrame(columns=
	['SAPAorthID','SAPAcov', 'SAPAnorm','SACEorthID','SACEcov', 'SACEnorm','GO','Note',
	'%Identity' ,'Gap_ID' ,'SAPAstart','SAPAend','SAPAchr','SACEstart','SACEend','SACEchr',
	'Ratio', 'Heterozygosity', 'Subtelomeric'])
	for i in range(len(sapaf)):
		sapaIter = list(sapaf.iterrows())[i][1]
		SAPAorthID = sapaIter.ID
		geneName = dictID[SAPAorthID]
		for gene in geneName:
			has_orth = SAPAorthID in orths.index
			if True:
				SACEorthID = gene
				SAPAcov = sapaIter.coverage
				try:
					SACEcov = sacef.loc[gene].coverage
				except:
					SACEcov = 0
				try:
					GO = sacedf.loc[gene].GO.replace(",",";")
				except:
					GO = "NA"
				try:
					Note = sacef.loc[gene].g_description.replace("%20"," ")
				except:
					Note = "NA"
				if SACEorthID != "NA":	#The next 2 lines are new, avoid miscalculationg Seq_ID
					SACEseq = "NA"		#when There is no SACEorthID
					if SACEorthID in seqIterSACE.keys():
						SACEseq = str(seqIterSACE[SACEorthID].seq)
					if SAPAorthID in seqIterSAPA.keys():
						SAPAseq = str(seqIterSAPA[SAPAorthID].seq)
					if SACEseq != "NA":
						res=list(align_sequences(SACEseq, SAPAseq))
						Seq_ID = res[1]
						Gap_ID = res[2]
					else:
						Seq_ID = 0
						Gap_ID = 0
				else:
					Seq_ID = 0
					Gap_ID = 0
				#This if was added temporaly (maybe) to study how we should define
				#orthologs #Said if was removed due to disuse 01/2024
				lssCoveredRds = 0
				mrCoveredRds = 0
				try:
					if SACEcov > SAPAcov:
						lssChr = sapadf.loc[SAPAorthID].seqname
						lssStrt = sapadf.loc[SAPAorthID].start
						lssEnd = sapadf.loc[SAPAorthID].end
						mrChr = sacedf.loc[gene].seqname
						mrStrt = sacedf.loc[gene].start
						mrEnd = sacedf.loc[gene].end
					else:
						lssChr = sacedf.loc[SAPAorthID].seqname
						lssStrt = sacedf.loc[SAPAorthID].start
						lssEnd = sacedf.loc[SAPAorthID].end
						mrChr = sapadf.loc[SAPAorthID].seqname
						mrStrt = sapadf.loc[SAPAorthID].start
						mrEnd = sapadf.loc[SAPAorthID].end
				except:
					if str(SACEcov) > (SAPAcov):
						lssChr = sapadf.loc[SAPAorthID].seqname
						lssStrt = sapadf.loc[SAPAorthID].start
						lssEnd = sapadf.loc[SAPAorthID].end
						mrChr = sacedf.loc[gene].seqname
						mrStrt = sacedf.loc[gene].start
						mrEnd = sacedf.loc[gene].end
					else:
						lssChr = sacedf.loc[SAPAorthID].seqname
						lssStrt = sacedf.loc[SAPAorthID].start
						lssEnd = sacedf.loc[SAPAorthID].end
						mrChr = sapadf.loc[SAPAorthID].seqname
						mrStrt = sapadf.loc[SAPAorthID].start
						mrEnd = sapadf.loc[SAPAorthID].end
				get_lss = f"samtools view {bam} '{lssChr}:{lssStrt}-{lssEnd}' | wc -l"
				get_mr = f"samtools view {bam} '{mrChr}:{mrStrt}-{mrEnd}' | wc -l"
				try:
					ps_lss = subprocess.Popen(get_lss,shell=True,stdout=subprocess.PIPE)
					ps_lss.wait()
					lessCoveredRds = float(str(ps_lss.communicate()[0]).split("'")[1].split("\\")[0])
				except:
					lessCoveredRds = 0
				try:
					ps_mr = subprocess.Popen(get_mr,shell=True,stdout=subprocess.PIPE)
					ps_mr.wait()
					mrCoveredRds = float(str(ps_mr.communicate()[0]).split("'")[1].split("\\")[0])
				except:
					mrCoveredRds = 0
				norm_file = f"/mnt/Timina/lmorales/Public/ymez/stats/coverage/{sample}/CONC/Q20/{sample}_CONC_Q20_perReference.csv"
				norm_open = pd.read_csv(norm_file, sep=",",names = ["Ref","Reads","Mean","Min","Max","l_quantile","Median","u_quantile","Above_10X"]).sort_values(by="Median",ascending=False)
				median = float(norm_open["Median"][1])
				SAPAnorm = float(SAPAcov)/median
				SACEnorm = float(SACEcov)/median
				ratio = SACEnorm/SAPAnorm
				if sapaIter.g_chr != 'SAPA_YPS138_v1_chr_mt':
					subt = subtelomeric(sapaIter.g_chr,sapaIter.g_start,sapaIter.g_end)
				else:
					subt = False
				if has_orth:
					if ratio >= 0.43 and ratio <= 2.3:
						het = "Het"
					elif ratio < 0.43:
						het = "HomSAPA"
					else:
						het = "HomSACE"
				else:
					if float(SAPAcov) > 0.85 * median:
						het = "HomSAPA"
					else:
						het = "HetSAPA"
				if lessCoveredRds > 0.25*mrCoveredRds and Seq_ID < 95:
					try:
						orthologs = orthologs.append({'SAPAorthID': SAPAorthID,
								'SAPAcov': SAPAcov, 'SAPAnorm': SAPAnorm,
								'SACEorthID': SACEorthID,
								'SACEcov': SACEcov, 'SACEnorm': SACEnorm,
								'GO': GO,'Note': Note,
								'%Identity' : Seq_ID,'Gap_ID' : Gap_ID,
								'SAPAstart':sapaIter.g_start,'SAPAend':sapaIter.g_end,'SAPAchr':sapaIter.g_chr,
								'SACEstart':sacef.loc[gene].g_start,'SACEend':sacef.loc[gene].g_end,
								'SACEchr':sacef.loc[gene].g_chr,
								'Ratio':ratio, 'Heterozygosity':het, 'Subtelomeric':subt,
								'Has_ortholog':has_orth},
								ignore_index=True)
					except:
						orthologs = orthologs.append({'SAPAorthID': SAPAorthID,
								'SAPAcov': SAPAcov, 'SAPAnorm': SAPAnorm,
								'SACEorthID': SACEorthID,
								'SACEcov': SACEcov, 'SACEnorm': SACEnorm,
								'GO': GO,'Note': Note,
								'%Identity' : Seq_ID,'Gap_ID' : Gap_ID,
								'SAPAstart':sapaIter.g_start,'SAPAend':sapaIter.g_end,'SAPAchr':sapaIter.g_chr,
								'SACEstart':"NA",'SACEend':"NA",
								'SACEchr':"NA",
								'Ratio':ratio, 'Heterozygosity':het, 'Subtelomeric':subt,
								'Has_ortholog':has_orth},
								ignore_index=True)
	return orthologs

def xrange(x):
    return iter(range(x))

def align_sequences(sequence_A, sequence_B, **kwargs):
    """
    Performs a global pairwise alignment between two sequences
    using the BLOSUM62 matrix and the Needleman-Wunsch algorithm
    as implemented in Biopython. Returns the alignment, the sequence
    identity and the residue mapping between both original sequences.
    """

    def _calculate_identity(sequenceA, sequenceB):
        """
        Returns the percentage of identical characters between two sequences.
        Assumes the sequences are aligned.
        """

        sa, sb, sl = sequenceA, sequenceB, len(sequenceA)
        matches = [sa[i] == sb[i] for i in xrange(sl)]
        seq_id = (100 * sum(matches)) / sl
        gapless_sl = sum([1 for i in xrange(sl) if (sa[i] != '-' and sb[i] != '-')])
        gap_id = (100 * sum(matches)) / gapless_sl
        return (seq_id, gap_id)
    matrix = kwargs.get('matrix', matlist.blosum62)
    gap_open = kwargs.get('gap_open', -10.0)
    gap_extend = kwargs.get('gap_extend', -0.5)
    alns = pairwise2.align.globalds(sequence_A, sequence_B,
                                    matrix, gap_open, gap_extend,
                                    penalize_end_gaps=(False, False) )
    best_aln = alns[0]
    aligned_A, aligned_B, score, begin, end = best_aln
    # Calculate sequence identity
    seq_id, g_seq_id = _calculate_identity(aligned_A, aligned_B)
    return ((aligned_A, aligned_B), seq_id, g_seq_id)
##
dictID = {}
sapadf =  pd.read_csv(sapatxt, sep="\t",names = ['seqname', 'start', 'end', 'attribute','ID','Notes'])
sacedf =  pd.read_csv(sacetxt, sep="\t",names = ['seqname', 'start', 'end', 'attribute','ID','Notes'])
orths = pd.read_csv(orths,sep="\t")
sacedf = sacedf.set_index("ID")
sapadf = sapadf.set_index("ID")
orths = orths.set_index("YPS138")
sacedf["GO"]  = sacedf.attribute.apply(lambda x: gene_info(x,"Ontology_term"))
for i,gene in sapadf.iterrows():
    dictID[gene_info(gene.attribute,"ID")] = gene_info(gene.attribute,"Name").split("/")
sapa = pd.DataFrame(csv2list(sapafile),columns = ["ID","g_chr","g_start","g_end","coverage","percAboveQMedian","g_description"])
sace = pd.DataFrame(csv2list(sacefile),columns = ["ID","g_chr","g_start","g_end","coverage","percAboveQMedian","g_description"])
sapa = sapa.sort_values(by = ["ID"]).reset_index().drop('index', axis=1)
sace = sace.sort_values(by = ["ID"]).reset_index().drop('index', axis=1).set_index("ID")
subt_r = pd.read_csv(subtelomeric_regions)
subt_r = subt_r.set_index("Chr")
o = find_orthologs(sapa,sace)
out_csv = out_folder +f"/{sample}_introSAPA_All_noOrth.csv"
o.to_csv(out_csv, index=False, header = True, sep=",")
