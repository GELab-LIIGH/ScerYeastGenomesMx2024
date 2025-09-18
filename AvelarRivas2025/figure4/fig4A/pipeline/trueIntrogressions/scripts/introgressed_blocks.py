import pandas as pd

file = snakemake.input[0]
annot = snakemake.input[1]
out = snakemake.output[0]


#This function groups the introgressions in blocks if they are consecutive
#introgressions.
def get_blocks(file,annot,out):
    #Open the introgressions file and the simplified annotation file from paradoxus.
    file_pd = pd.read_csv(file).sort_values(by=["SAPAchr","SAPAstart"]).drop_duplicates(subset="SAPAorthID").reset_index(drop=True)
    annot_pd = pd.read_csv(annot,sep="\t",header=None).set_axis(["seqname","start","end","attribute","ID","comments"],axis=1).sort_values(by=["seqname","start"]).reset_index(drop=True)
    #Define necessary empty lists which will be used to create the dataframe with
    #the blocks
    size_l = []
    bname_l = []
    sizeNt_l = []
    chr_l = []
    start_l = []
    end_l = []
    genes_l = []
    #Counter which tells in what introgressed gene we are currently at.
    current_idx = 0
    #To avoid an out-of-bounds error current_idx goes from 0 to len(file_pd)-1.
    while current_idx < file_pd.index.stop-1:
        #Get id and chromosome from the introgression that will start the block,
        #the block seed, which can be the only member of the block.
        #In the begining th block seed and the last member of the block are the
        #same.
        current_gene = file_pd.iloc[current_idx]["SAPAorthID"]
        current_chr = file_pd.iloc[current_idx]["SAPAchr"]
        #Size of the current block, whose seed is current_gene.
        #Block sizes start at 1 since a block originally is comprised by only one
        #gene, the current_gene, it is possible for it to stay at size one.
        size = 1
        #next_idx is the index used to see if the introgression next to the current
        #block is the next gene in the annotation file of paradoxus.
        next_idx = current_idx+size
        #Get the id and chromosome of the introgression next to the current block,
        #which is candidate to also be part of such block.
        next_gene = file_pd.iloc[next_idx]["SAPAorthID"]
        next_chr = file_pd.iloc[next_idx]["SAPAchr"]
        #annot_idx is the index from the annot_pd dataframe row which corresponds
        #to the block seed.
        annot_idx = annot_pd.index[annot_pd["ID"]==current_gene].values[0]
        #start and end of the block will be extracted from the annotation file.
        #start is the start coords of the block seed.
        start = annot_pd.iloc[annot_idx]["start"]
        #end is the end of the last member of the block, originally it is the end
        #of the seed of the block.
        end = annot_pd.iloc[annot_idx]["end"]
        #sizeNT is the size of the block in nucleotides, originally the sizeNT is
        #the same size as the gene introgressed
        sizeNT = end - start
        #genes contain the genes that conform the block, separated by semicolon,
        #i.e.:gene1;gene2;gene3. Originally it is equal to the current_gene
        genes = current_gene
        #next is the index of the dataframe row that contains the gene bellow the
        #last member of the current block, which is candidate to also be part of
        #the block.
        next = annot_idx+size
        #to_match is the gene bellow the last member in the current block IN THE
        #PARADOXUS ANNOTATION, and is candidate to be the next member of the current
        #block, thus the name. Originally to_match is the gene immediatly bellow
        #annot_pd.iloc[annot.idx]
        to_match = annot_pd.iloc[next]["ID"]
        #match is a boolean, false if to_match is different to next_gene or
        #current_chr is different from next_chr.
        match = to_match == next_gene and current_chr == next_chr
        #If match is false, the block is broken and the information of the current
        #block will be added to the lists defined above (while wont be executed).
        #If match is true, the while wil be executing, updating the information
        #of the block, the blocks size and information will be updated until the
        #match equals false.
        while match:
            #block size is incremented by one since match == true, meaning the
            #block contains one more gene now.
            size += 1
            #End of block is updated with the end of the new last member of the
            #current block.
            end = annot_pd.iloc[next]["end"]
            #Once the end is update, the size of the block in nucleotides is also
            #updated.
            sizeNT = end - start
            #The new last member of the block is added to the genes list.
            genes = genes + ";"+ to_match
            #The information of the new introgression cadidate to be member of the
            #current block is updated.
            next_idx = current_idx+size
            try:
                next_gene = file_pd.iloc[next_idx]["SAPAorthID"]
            except IndexError:
                next_gene = "ended"
            try:
                next_chr = file_pd.iloc[next_idx]["SAPAchr"]
            except IndexError:
                next_chr = "ended"
            #Information related to the annot_pd is updated
            next = annot_idx+size
            to_match = annot_pd.iloc[next]["ID"]
            #match is updated
            match = to_match == next_gene and current_chr == next_chr
        bname = f"{genes.split(';')[0]}_2_{genes.split(';')[-1]}"
        #The block is now broken, updating the lists defined above.
        size_l.append(size)
        bname_l.append(bname)
        sizeNt_l.append(sizeNT)
        chr_l.append(current_chr)
        start_l.append(start)
        end_l.append(end)
        genes_l.append(genes)
        #current_idx is updated
        current_idx +=size
    #Since the while avoids entering the last row of the file_pd dataframe, if the
    #last introgression isn't a member of the last block, the last gene wont be
    #included in the bloks dataframe, as its own block of size one. Thus the
    #following block of code.
    if not match and next_gene != "ended" :
        current_gene = file_pd.iloc[current_idx]["SAPAorthID"]
        annot_idx = annot_pd.index[annot_pd["ID"]==current_gene].values[0]
        size_l.append(1)
        bname_l.append(f"{current_gene}_2_{current_gene}")
        start = annot_pd.iloc[annot_idx]["start"]
        end = annot_pd.iloc[annot_idx]["end"]
        sizeNt = end - start
        sizeNt_l.append(sizeNt)
        chr_l.append(file_pd.iloc[current_idx]["SAPAchr"])
        start_l.append(start)
        end_l.append(end)
        genes_l.append(current_gene)
    #The blocks dataframe is created with the lists defined above
    blocks = pd.DataFrame({"Name":bname_l,"Size":size_l,"SizeNT":sizeNt_l,"Chr":chr_l,"Start":start_l,"End":end_l,"Genes":genes_l})
    #The dataframe is written as a csv
    blocks.to_csv(out, index=False, header = True, sep=",")

    return blocks



if __name__ == '__main__':
    get_blocks(file,annot,out)
