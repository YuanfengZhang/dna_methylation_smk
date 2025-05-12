configfile: "config/runtime_config.yaml"
from textwrap import dedent


# ! not working, See https://github.com/smithlabcode/dnmtools/issues/260

rule dnmtools_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        temp("result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.dnmtools.tmp"),
        "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.dnmtools.gz"
    benchmark:
        "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.count.benchmark"
    params:
        ref            = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        extra_params   = config["dnmtools"]["count"]["extra_params"] or ""
    threads: 16
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.CountParentDir}
        mkdir -p dnmtools
        if [ ! -f {wildcards.BaseName}.sorted.bam ]; then
            samtools sort -@ {threads} -o {wildcards.BaseName}.sorted.bam {wildcards.BaseName}.bam
            samtools index -@ {threads} {wildcards.BaseName}.sorted.bam
        fi
        dnmtools counts \\
            -c {params.ref} \\
            -t {threads} {params.extra_params} \\
            -o dnmtools/{wildcards.BaseName}.dnmtools.tmp \\
            {wildcards.BaseName}.sorted.bam
        cd dnmtools

        awk -F'\\t' '$6 > 4' \\
            {wildcards.BaseName}.dnmtools.tmp \\
            | pigz  -p {threads} --best > {wildcards.BaseName}.dnmtools.gz
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
