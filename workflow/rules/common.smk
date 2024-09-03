
import pandas as pd
import os

configfile: "config/config.yaml"

ruleorder: generate_candidates > preprocess

# construct genome name
datatype_genome = "dna"
species = config["ref"]["species"]
build = config["ref"]["build"]
release = config["ref"]["release"]
genome_name = f"genome.{datatype_genome}.{species}.{build}.{release}"
genome_prefix = f"resources/{genome_name}"
genome = f"{genome_prefix}.fasta"
genome_fai = f"{genome}.fai"

# list of alleles that will be checked
loci = ["A","B", "C", "DQA1","DQB1"]

# read samples
samples = pd.read_csv(config["samples"], sep="\t").set_index(
    ["sample_name"], drop=False
)

def get_fastq_input(wildcards):
    sample = samples.loc[wildcards.sample]
    return [sample["fq1"], sample["fq2"]]

