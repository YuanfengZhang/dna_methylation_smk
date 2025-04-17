configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule astair_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        multiext("result/{BaseName}/{CountParentDir}"
                 "/astair/{BaseName}",
                 ".astair.mods.gz",
                 ".astair.stats")
    benchmark:
        "result/{BaseName}/{CountParentDir}/astair/{BaseName}.count.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["astair"][wildcards.BaseName.split('_')[1]],
        method       = lambda wildcards: config["astair"]["method"][wildcards.BaseName.split('_')[0][: 2]],
        extra_params = config["astair"]["count"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.CountParentDir}
        astair call \
            -i {wildcards.BaseName}.bam \
            --reference {params.ref} \
            -co CpG -ni . -d astair -t {threads} \
            --method {params.method} {params.extra_params}
        cd astair
        mv \
            {wildcards.BaseName}_{params.method}_CpG.mods \
            {wildcards.BaseName}.astair.mods
        pigz --best {wildcards.BaseName}.astair.mods
        mv \
            {wildcards.BaseName}_{params.method}_CpG.stats \
            {wildcards.BaseName}.astair.stats
        """)
