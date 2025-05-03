configfile: "config/runtime_config.yaml"
from textwrap import dedent

"""
Although mergeContext only utilize 1 CPU, it has to be merged in a single rule with extract,
otherwise the snakemake will weirdly check the deleted input and throw an error.
I don't want talk too much about it. The extract always finishes in a few minutes,
so it won't waste too much resources.
"""


rule methyldackel_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        "result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.bedgraph.gz",
        "result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.merged.bedgraph.gz"
    benchmark:
        "result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.count.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        mode         = lambda wildcards: "-mCtoT" if wildcards.BaseName.split('_')[0][: 2] == "PS" else "",
        extra_params = config["methyldackel"]["count"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.CountParentDir}
        mkdir -p methyldackel
        MethylDackel extract \\
            {params.ref} {wildcards.BaseName}.bam \\
            -o methyldackel/{wildcards.BaseName} \\
            -@ {threads} {params.mode} {params.extra_params}

        cd methyldackel

        mv \\
            {wildcards.BaseName}_CpG.bedGraph \\
            {wildcards.BaseName}.bedgraph
        
        MethylDackel mergeContext \\
            {params.ref} {wildcards.BaseName}.bedgraph \\
            -o {wildcards.BaseName}.merged.bedgraph {params.extra_params}
        
        sed -i '1d' *.bedgraph

        pigz --best -p 8 {wildcards.BaseName}.bedgraph
        pigz --best -p 8 {wildcards.BaseName}.merged.bedgraph
        """)
