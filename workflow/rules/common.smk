
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
samples = pd.read_csv(config["orthanq_input"], sep="\t").set_index(
    ["sample_name"], drop=False
)

#The workflow will by default use BAM input if 'bam' column is filled. Otherwise, it will use `fq1` and `fq2` columns.
def get_fastq_input(wildcards):
    sample = samples.loc[wildcards.sample]
    return [sample["fq1"], sample["fq2"]]

def get_orthanq_input(wildcards):
    sample = samples.loc[wildcards.sample]
    if pd.notna(sample.get("bam")) and sample["bam"]:
        return sample["bam"]
    return get_fastq_input(wildcards)

def get_orthanq_input_params(wildcards):
    sample = samples.loc[wildcards.sample]
    if pd.notna(sample.get("bam")) and sample["bam"]:
        return f"--bam-input {sample['bam']}"
    fq1, fq2 = get_fastq_input(wildcards)
    return f"--reads {fq1} {fq2}"