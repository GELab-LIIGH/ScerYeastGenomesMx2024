import pandas as pd

samples = snakemake.input
out = snakemake.output[0]

def get_strain(path):
    strain = path.split("/")[-1].split("_")[0]
    return strain

strains_l = []
no_l = []
intros_l = []
avg_l = []
sd_l = []
min_l = []
max_l = []
q25_l = []
q50_l = []
q75_l = []
avgNT_l = []
sdNT_l = []
minNT_l = []
maxNT_l = []
q25NT_l = []
q50NT_l = []
q75NT_l = []

for s in samples:
    s_pd = pd.read_csv(s)
    strains_l.append(get_strain(s))
    intros_l.append(s_pd["Size"].sum())
    no_l.append(len(s_pd))
    avg_l.append(s_pd["Size"].mean())
    sd_l.append(s_pd["Size"].std())
    max_l.append(s_pd["Size"].max())
    min_l.append(s_pd["Size"].min())
    q25_l.append(s_pd["Size"].quantile(0.25))
    q50_l.append(s_pd["Size"].quantile(0.50))
    q75_l.append(s_pd["Size"].quantile(0.75))
    avgNT_l.append(s_pd["SizeNT"].mean())
    sdNT_l.append(s_pd["SizeNT"].std())
    maxNT_l.append(s_pd["SizeNT"].max())
    minNT_l.append(s_pd["SizeNT"].min())
    q25NT_l.append(s_pd["SizeNT"].quantile(0.25))
    q50NT_l.append(s_pd["SizeNT"].quantile(0.50))
    q75NT_l.append(s_pd["SizeNT"].quantile(0.75))


summary = pd.DataFrame({"Strain":strains_l,"NoIntrogressions":intros_l,
    "No_Blocks":no_l,"Avg_Size":avg_l,"Sd_Size":sd_l,"Max_Size":max_l,
    "Min_Size":min_l,"Q25_Size":q25_l,"Q50_Size":q50_l,"Q75_Size":q75_l,
    "No_Blocks":no_l,"Avg_SizeNT":avgNT_l,"Sd_SizeNT":sdNT_l,"Max_SizeNT":maxNT_l,
    "Min_SizeNT":minNT_l,"Q25_SizeNT":q25NT_l,"Q50_SizeNT":q50NT_l,"Q75_SizeNT":q75NT_l})
summary.to_csv(out, index=False, header = True, sep=",")
