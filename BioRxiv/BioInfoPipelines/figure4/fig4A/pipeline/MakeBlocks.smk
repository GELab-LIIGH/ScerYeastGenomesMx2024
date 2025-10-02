#Configuration file with defined variables
configfile: "config.yaml"

#Parse arguments
introgressions = config["introgressions_folder"]
annotation = config["annotation_file"]
subsample = config["subsample"]
samples = config["sace"]

rule make_blocks:
    input:
        file = introgressions + "{sample}_introSAPA_All_noOrth.csv",
        annot = annotation
    output:
        out = "data/blocks/{sample}_blocks.csv"
    script:
        "scripts/introgressed_blocks.py"

rule get_statistics:
    input:
        strains = expand("data/blocks/{sample}_blocks.csv",sample=samples)
    output:
        "data/summary_statistics_blocks.csv"
    script:
        "scripts/block_summary.py"
