configfile: "config/runtime_config.yaml"

rule methyldackel_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        "result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.bedgraph"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        mode         = lambda wildcards: "-mCtoT" if wildcards.BaseName.split('_')[0] == "PS" else "",
        extra_params = (config["methyldackel"]["count"]["extra_params"]
                        if config["methyldackel"]["count"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.BaseName}/{wildcards.CountParentDir}
        MethylDackel extract \
            {params.ref} {wildcards.BaseName}.bam \
            -o methyldackel/{wildcards.BaseName} \
            -@ {threads} {params.mode} {params.extra_params}
        mv \
            methyldackel/{wildcards.BaseName}_CpG.bedGraph \
            methyldackel/{wildcards.BaseName}.bedgraph
        """

rule methyldackel_merge_context:
    input:
        ("result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.bedgraph")
    output:
        "result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.bedgraph.gz",
        "result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.merged.bedgraph.gz"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        extra_params = (config["methyldackel"]["count"]["extra_params"]
                        if config["methyldackel"]["count"]["extra_params"] else "")
    threads: 1
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.BaseName}/{wildcards.CountParentDir}/methyldackel
        MethylDackel mergeContext \
            {params.ref} {wildcards.BaseName}.bedgraph \
            -o {wildcards.BaseName}.merged.bedgraph {params.extra_params}
        
        sed -i '1d' *.bedgraph

        pigz --best -p 8 {wildcards.BaseName}.bedgraph
        pigz --best -p 8 {wildcards.BaseName}.merged.bedgraph
        """
