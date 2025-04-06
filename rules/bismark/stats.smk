configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule bam2nuc:
    input:
        "result/{BaseName}/{StatsParentDir}/{BaseName}.bam"
    output:
        "result/{BaseName}/{StatsParentDir}/{BaseName}.bam2nuc.txt"
    params:
        ref          = lambda wildcards: config["ref"]["bismark-bowtie2"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bismark"]["bam2nuc"]["extra_params"] or ""
    threads: 1
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.StatsParentDir}
        bam2nuc \
            {wildcards.BaseName}.bam \
            --dir . \
            --genome_folder {params.ref} \
            {params.extra_params}
        mv \
            {wildcards.BaseName}.nucleotide_stats.txt \
            {wildcards.BaseName}.bam2nuc.txt
        """)
