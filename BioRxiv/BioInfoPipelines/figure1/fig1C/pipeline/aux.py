"""
Auxiliary functions for the snakemake pipeline
"""
import pandas as pd
from glob import glob


def get_samples():
	files = glob("data/*_CONC.g.vcf")
	samples = [s.split("/")[1].split("_")[0] for s in files]
	return samples
