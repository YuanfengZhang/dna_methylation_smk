configfile: "config/runtime_config.yaml"

rule methyldackel_mbias:
    input:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}bam"
    output:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}methydackel_mbias_OT.svg",
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}methydackel_mbias_OB.svg"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        extra_params = (config["methyldackel"]["mbias"]["extra_params"]
                        if config["methyldackel"]["mbias"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.BaseName}/{wildcards.BamStatsParentDir}
        MethylDackel mbias \
            {params.ref} {wildcards.BaseName}.bam \
            ./{wildcards.BaseName}.methydackel_mbias \
            -@ {threads} {params.extra_params}
        """
