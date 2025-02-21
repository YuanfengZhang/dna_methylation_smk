configfile: "config/runtime_config.yaml"

rule methyldackel_count:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        ("result/{fname}/{trimmer}/{aligner}/{deduper}/"
         "methyldackel/{fname}{bqsr}bedgraph")
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.fname.split('_')[1]],
        extra_params = (config["methyldackel"]["count"]["extra_params"]
                        if config["methyldackel"]["count"]["extra_params"] else ""),
        pattern      = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        MethylDackel extract \
            {params.ref} {wildcards.fname}{params.pattern}.bam \
            -o methyldackel/{wildcards.fname}{params.pattern} \
            -@ {threads} {params.extra_params}
        mv \
            methyldackel/{wildcards.fname}{params.pattern}_CpG.bedGraph \
            methyldackel/{wildcards.fname}{params.pattern}.bedgraph
        """

rule methyldackel_merge_context:
    input:
        ("result/{fname}/{trimmer}/{aligner}/{deduper}/methyldackel/{fname}{bqsr}bedgraph")
    output:
        ("result/{fname}/{trimmer}/{aligner}/{deduper}/"
         "methyldackel/{fname}{bqsr}bedgraph.gz"),
        ("result/{fname}/{trimmer}/{aligner}/{deduper}/"
         "methyldackel/{fname}{bqsr}merged.bedgraph.gz")
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.fname.split('_')[1]],
        extra_params = (config["methyldackel"]["count"]["extra_params"]
                        if config["methyldackel"]["count"]["extra_params"] else ""),
        pattern      = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 1
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}/methyldackel
        MethylDackel mergeContext \
            {params.ref} {wildcards.fname}{params.pattern}.bedgraph \
            -o {wildcards.fname}{params.pattern}.merged.bedgraph {params.extra_params}

        pigz --best {wildcards.fname}{params.pattern}.bedgraph
        pigz --best {wildcards.fname}{params.pattern}.merged.bedgraph
        """
