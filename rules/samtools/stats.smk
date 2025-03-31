configfile: "config/runtime_config.yaml"

rule samtools_stats:
    input:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}bam"
    output:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}samtools.stats.txt"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        extra_params = (config["samtools"]["stats"]["extra_params"]
                        if config["samtools"]["stats"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        samtools stats {input} \
            -d \
            -r {params.ref} \
            -@ {threads} \
            {params.extra_params} \
            > {output}
        """

rule samtools_flagstat:
    input:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}bam"
    output:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}samtools.flagstats.txt"
    params:
        extra_params = (config["samtools"]["flagstat"]["extra_params"]
                        if config["samtools"]["flagstat"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        samtools flagstat \
            -@ {threads} \
            -O tsv \
            {input} > {output}
        """
