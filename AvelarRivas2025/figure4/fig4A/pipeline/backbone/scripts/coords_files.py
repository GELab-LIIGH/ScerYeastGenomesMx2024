import pandas as pd
#Strains, passed as list from Snakemake config file
samples = snakemake.config["sace"]
folder = snakemake.config["blocks_folder"]
out_coords = snakemake.output[0]
out_shared = snakemake.output[1]
#Function that gets the file with the introgressions of each
#strain and returns a dict with strain:strain_introgression_file
#pairs
def paths(strains):
    paths_dict = {}
    for s in strains:
        paths_dict[s] = folder + s +"_blocks.csv"
    return paths_dict

def update_universe(intros_df,strain,universe,shared):
    name = intros_df["Name"]
    chromosome = intros_df["Chr"]
    start = intros_df["Start"]
    end = intros_df["End"]
    if name not in universe.keys():
        universe[name] = [chromosome,start,end]
    else:
        pass
    if name not in shared.keys():
        shared[name] = [strain]
    else:
        if strain not in shared[name]: ##Por si acaso
            shared[name].append(strain)
    return name

def get_intros(strains,folder,out1,out2):
    intro_universe = dict()
    files = paths(strains)
    shared_dict = dict()
    for strain in strains:
        print(strain)
        reader = pd.read_csv(files[strain])
        for index, intro in reader.iterrows():
            update_universe(intro,strain,intro_universe,shared_dict)
        del reader
    print(len(intro_universe.keys()))
    df_coords = pd.DataFrame.from_dict(intro_universe,orient="index",columns=["Chromosome","Start","End"])
    df_coords.index.name = "ID"
    df_coords.to_csv(out1,index = True,header = True,sep = "\t")
    df_shared = pd.DataFrame({'Gene': list(shared_dict.keys()), 'Id': list(shared_dict.values())})
    df_shared = df_shared.explode(column='Id').reset_index(drop=True)
    df_shared.to_csv(out2,index=False,header=True,sep="\t")
    return (out1,out2)

if __name__ == '__main__':
    get_intros(samples,folder,out_coords,out_shared)
