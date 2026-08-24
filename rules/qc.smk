rule sylph:
    input:
       fastqs=expand(["data/trimmed/{sample}_1.fastq.gz", "data/trimmed/{sample}_2.fastq.gz"], sample=samples["sample"]),
       db="data/syl_db/gtdb-r220-c200-dbv1.syldb"
    output:
        "data/qc/sylph_profiles.tsv"
    conda:
        "../envs/sylph.yml"
    threads: 12
    resources:
        runtime=lambda wildcards, input: 5 + (len(input) * 0.01),
        mem_mb=32000
    shell:
        "sylph profile {input.db} -1 data/trimmed/*_1.fastq.gz -2 data/trimmed/*_2.fastq.gz -t {threads} > {output}"

rule filter_tophits:
    input:
        "data/qc/sylph_profiles.tsv"
    output:
        "data/qc/sylph_profiles_tophits.tsv"
    shell:
        """
        python3 scripts/filter_tophits.py {input} {output}
        """

rule summarize_fastp:
    input:
        expand("data/qc/fastp/{sample}.json", sample=samples["sample"])
    output:
        "data/qc/fastp_summary.txt"
    shell:
        """
        python scripts/summarize_fastp.py data/qc/fastp
        """

checkpoint filter_read_qc_metrics:
    input:
        sylph="data/qc/sylph_profiles_tophits.tsv",
        fastp="data/qc/fastp_summary.txt"
    params:
        genus=config["genus"],
        species=config["species"],
        min_coverage=config["min_q30_coverage"],
        genome_size=config["expected_genome_size"],
        min_abundance=config["min_tax_read_assignment"]
    output:
        combined_metrics="data/qc/read_metrics.tsv",
        pass_samples="data/qc/read_metrics_pass_samples.txt"
    shell:
        """
        python scripts/apply_read_qc_thresholds.py --genus {params.genus} --species {params.species} --genome_size {params.genome_size} --min_coverage {params.min_coverage} --min_abundance {params.min_abundance} {input.sylph} {input.fastp} {output.combined_metrics} {output.pass_samples}
        """

rule create_checkm2_input:
    input:
        "data/annotations/{sample}/{sample}.faa"
    output:
        temp("data/qc/checkm2_input/{sample}.faa")
    shell:
        """
        mkdir -p data/qc/checkm2_input/
        ln -rs {input} {output}
        """

def get_checkm2_input(wc):
    return expand("data/qc/checkm2_input/{sample}.faa", sample=get_passed_reads(wc))

rule checkm2:
    input:
        proteins=get_checkm2_input,
        db="data/CheckM2_database/uniref100.KO.1.dmnd"
    output:
        report="data/qc/checkm2/quality_report.tsv",
    threads: 16
    conda:
        "../envs/checkm2.yml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 1440,
        mem_mb=32000
    shell:
        """
        checkm2 predict --force -t {threads} -x .faa --genes --output-directory data/qc/checkm2/ --input data/qc/checkm2_input/ --database_path {input.db}
        """

def get_quast_input(wc):
    return expand("data/filtered_assemblies/{sample}_contigs_filtered.fa", sample=get_passed_reads(wc))

rule quast:
    input:
        assemblies=get_quast_input
    output:
        "data/qc/quast/transposed_report.tsv"
    threads: 16
    conda:
        "../envs/quast.yml"
    resources:
        runtime=lambda wildcards, attempt: attempt * 30,
        mem_mb=8000
    shell:
        """
        quast -o data/qc/quast --threads 16 --no-plots {input.assemblies}
        """
