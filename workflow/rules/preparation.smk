rule get_genome:
    output:
        genome,
    log:
        "logs/preparation/get-genome.log",
    params:
        species=config["ref"]["species"],
        datatype="dna",
        build=config["ref"]["build"],
        release=config["ref"]["release"],
        chromosome=config["ref"].get("chromosome"),
    cache: "omit-software"
    wrapper:
        "v2.3.2/bio/reference/ensembl-sequence"


rule genome_faidx:
    input:
        genome,
    output:
        genome_fai,
    log:
        "logs/preparation/genome-faidx.log",
    cache: "omit-software"
    wrapper:
        "v2.3.2/bio/samtools/faidx"


rule get_hla_genes_and_xml:
    output:
        genes="results/preparation/hla_gen.fasta",
        xml="results/preparation/hla.xml.zip",
    log:
        "logs/get_hla_genes_and_xml.log",
    params:
        genes_link="ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta",
        xml_link="ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/xml/hla.xml.zip",
    shell:
        "wget -c  {params.genes_link} -O {output.genes} && "
        "wget -c {params.xml_link} -O {output.xml} 2> {log}"


rule unzip_xml:
    input:
        "results/preparation/hla.xml.zip",
    output:
        xml="results/preparation/hla.xml",
    log:
        "logs/unzip_xml.log",
    params:
        path_to_unzip=lambda wc, output: os.path.dirname(output.xml),
    shell:
        "unzip -o {input} -d {params.path_to_unzip}"


rule get_allele_frequencies:
    output:
        rda="results/preparation/allele_frequencies.rda",
        csv="results/preparation/allele_frequencies.csv",
    conda:
        "../envs/R.yaml"
    log:
        "logs/get_allele_frequencies.log",
    params:
        data_link="https://github.com/Genentech/midasHLA/blob/11bde30cbbf11b34f2dea29a6284371a9c1e9440/data/allele_frequencies.rda?raw=true",
    shell:
        """
        Rscript -e 'download.file("{params.data_link}", "{output.rda}"); load("{output.rda}"); write.csv(allele_frequencies, file="{output.csv}")' 2> {log}
        """


rule get_pangenome:
    output:
        "results/preparation/hprc-v1.0-mc-grch38.xg",
    log:
        "logs/get_pangenome.log",
    params:
        data_link="https://s3-us-west-2.amazonaws.com/human-pangenomics/pangenomes/freeze/freeze1/minigraph-cactus/hprc-v1.0-mc-grch38.xg",
    shell:
        "wget -c {params.data_link} -O {output} 2> {log}"


rule bwa_index:
    input:
        genome,
    output:
        idx=multiext(genome, ".amb", ".ann", ".bwt", ".pac", ".sa"),
    log:
        "logs/preparation/bwa_index.log",
    cache: True
    wrapper:
        "v2.3.2/bio/bwa/index"
