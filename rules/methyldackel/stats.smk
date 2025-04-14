configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule methyldackel_mbias:
    input:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.bam"
    output:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.methydackel_mbias_OT.svg",
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.methydackel_mbias_OB.svg"
    benchmark:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.methydackel_mbias.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        extra_params = config["methyldackel"]["mbias"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.BamStatsParentDir}
        MethylDackel mbias \
            {params.ref} {wildcards.BaseName}.bam \
            ./{wildcards.BaseName}.methydackel_mbias \
            -@ {threads} {params.extra_params}
        """)
