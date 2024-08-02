
import pandas as pd
import os

configfile: "config/config.yaml"

ruleorder: generate_candidates > preprocess

# list of alleles that will be checked
loci = ["A","B","C","DQA1","DQB1"]

# read samples
samples = pd.read_csv(config["samples"], sep="\t").set_index(
    ["sample_name"], drop=False
)

def get_fastq_input(wildcards):
    sample = samples.loc[wildcards.sample]
    return [sample["fq1"], sample["fq2"]]