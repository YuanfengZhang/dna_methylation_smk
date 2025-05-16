configfile: "config/runtime_config.yaml"
from textwrap import dedent


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
