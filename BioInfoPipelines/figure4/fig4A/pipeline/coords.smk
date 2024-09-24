#Configuration file with defined variables
configfile: "config.yaml"
#Parse arguments
coords = config["coords"]
shared = config["shared"]

rule coords_files:
    output:
        coords,
        shared
    script:
        "scripts/coords_files.py"
