configfile: "config/runtime_config.yaml"
from textwrap import dedent


def get_method_param(wildcards):
    match wildcards.BaseName.split('_')[0][: 2]:
    case "PS":
        return ""
    case _:
        return "--method CtoT"


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
        method       = lambda wildcards: get_method_param(wildcards),
        extra_params = config["astair"]["count"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.CountParentDir}
        astair call \\
            -i {wildcards.BaseName}.bam \\
            --reference {params.ref} \\
            --context CpG -ni . -d astair -t {threads} \\
            --method {params.method} {params.extra_params}
        cd astair
        mv \\
            {wildcards.BaseName}_{params.method}_CpG.mods \\
            {wildcards.BaseName}.astair.mods
        pigz --best {wildcards.BaseName}.astair.mods
        mv \\
            {wildcards.BaseName}_{params.method}_CpG.stats \\
            {wildcards.BaseName}.astair.stats
        """)
