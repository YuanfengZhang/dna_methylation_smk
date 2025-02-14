configfile: "config/runtime_config.yaml"

rule methyldackel_mbias:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        ("result/{fname}/{trimmer}/{aligner}/{deduper}/"
         "{fname}{bqsr}methydackel_mbias_OT.svg")
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.fname.split('_')[1]],
        extra_params = (config["methyldackel"]["mbias"]["extra_params"]
                        if config["methyldackel"]["mbias"]["extra_params"] else ""),
        pattern      = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        MethylDackel mbias \
            {params.ref} {wildcards.fname}{params.pattern}.bam \
            ./{wildcards.fname}{params.pattern}.methydackel_mbias {params.extra_params}
        """
