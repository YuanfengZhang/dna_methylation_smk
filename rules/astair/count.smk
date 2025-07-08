configfile: "config/runtime_config.yaml"
from textwrap import dedent


def get_method_param(wildcards):
    match wildcards.BaseName.split('_')[0][: 2]:
        case "PS":
            return ""
        case _:
            return "--method CtoT"


def get_method_param2(wildcards):
    match wildcards.BaseName.split('_')[0][: 2]:
        case "PS":
            return "_mCtoT"
        case _:
            return "_CtoT"


rule astair_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        mods         = "result/{BaseName}/{CountParentDir}/astair/{BaseName}.astair.mods.gz",
        stats        = "result/{BaseName}/{CountParentDir}/astair/{BaseName}.astair.stats"
    benchmark:
        "result/{BaseName}/{CountParentDir}/astair/{BaseName}.count.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        method       = lambda wildcards: get_method_param(wildcards),
        method2      = lambda wildcards: get_method_param2(wildcards),
        extra_params = config["astair"]["count"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        export tmp_dir="result/{wildcards.BaseName}/{wildcards.CountParentDir}"
        mkdir -p "${{tmp_dir}}/astair"
        astair call \\
            -i "${{tmp_dir}}/{wildcards.BaseName}.bam" \\
            --reference {params.ref} \\
            --context CpG -ni . -d "${{tmp_dir}}/astair" \\
            -t {threads} {params.method} {params.extra_params}

        pigz --best -p {threads} -c \\
          "${{tmp_dir}}/astair/{wildcards.BaseName}{params.method2}_CpG.mods" \\
          > {output.mods}
        mv \\
            "${{tmp_dir}}/astair/{wildcards.BaseName}{params.method2}_CpG.stats" \\
            {output.stats}
        rm -f "${{tmp_dir}}/astair/{wildcards.BaseName}{params.method2}_CpG.mods"
        """)
