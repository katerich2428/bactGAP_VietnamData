import glob
import pandas as pd
from Bio import SeqIO

configfile: "config/config.yaml"

samples = pd.read_csv(config["samples"], sep="\t")
samples_dict = samples.set_index("sample").to_dict(orient="index")

def get_passed_reads(wc):
    checkpoint_output = checkpoints.filter_read_qc_metrics.get(**wc).output.pass_samples
    df = pd.read_csv(checkpoint_output, sep="\t")
    return df["sample"].tolist()

include: "rules/download.smk"
#include: "rules/reference_mapping.smk"
#include: "rules/trees.smk"
include: "rules/assembly_annotation.smk"
include: "rules/amr.smk"
#include: "rules/clustering.smk"
include: "rules/qc.smk"
#include: "rules/salmonella.smk"

localrules: all, filter_tophits, filter_read_qc_metrics, create_checkm2_input, summarize_fastp

rule all:
    input:
        #expand("data/annotations/{sample}/{sample}.gff3", sample=samples["sample"]),
        "data/qc/checkm2/quality_report.tsv",
        "data/qc/quast/transposed_report.tsv",
        expand("data/amrfinder/{sample}.txt", sample=samples["sample"]),
        #"data/alignments/core.full.aln",
        #"data/trees/gubbins/core.final_tree.tre",
        #"data/poppunk/poppunk_clusters.csv",
        #expand("data/seqsero2/{sample}/SeqSero_result.tsv", sample=samples["sample"]) if config["genus"] == "Salmonella" else [],





