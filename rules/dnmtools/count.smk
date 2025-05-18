configfile: "config/runtime_config.yaml"
from textwrap import dedent


# ! the latest released dnmtools cannot skip the reads with $3 (chromosome) that not in the ref.fa,
# ! See https://github.com/smithlabcode/dnmtools/issues/260
# ! The current solution is to use awk to get the legal records, run dnmtools count, and then remove it.

rule dnmtools_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        temp_sam          = temp("result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.mapped.sam"),
        temp_bedgraph     = temp("result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.dnmtools.tmp"),
        bedgraph          = "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.dnmtools.gz"
    benchmark:
        "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.count.benchmark"
    params:
        ref               = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        extra_params      = config["dnmtools"]["count"]["extra_params"] or ""
    threads: 16
    conda:
        "conda.yaml"
    shell:
        dedent("""
        dnmtools_dir="result/{wildcards.BaseName}/{wildcards.CountParentDir}/dnmtools"
        mkdir -p "$dnmtools_dir"
        if [ ! -f {output.temp_sam} ]; then
            samtools view -h {input} |\\
            awk '$3 != "*"' > {output.temp_sam}
        fi
        dnmtools counts \\
            -c {params.ref} -t {threads} \\
            -progress {params.extra_params} \\
            -o {output.temp_bedgraph} \\
            {output.temp_sam}

        awk -F'\\t' '$6 > 4' \\
            {output.temp_bedgraph} \\
            | pigz -p {threads} --best > {output.bedgraph}
        """)


rule dnmtools_levels:
    input:
        "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.dnmtools.tmp"
    output:
        "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.dnmtools.stats"
    benchmark:
        "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.level.benchmark"
    params:
        extra_params   = config["dnmtools"]["levels"]["extra_params"] or ""
    threads: 16
    conda:
        "conda.yaml"
    shell:
        dedent("""
        dnmtools levels -o {output} {input} {params.extra_params}
        """)
