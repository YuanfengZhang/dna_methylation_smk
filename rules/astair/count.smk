configfile: "config/runtime_config.yaml"

rule astair_count:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        multiext("result/{fname}/{trimmer}/{aligner}/{deduper}"
                 "/astair/{fname}{bqsr}",
                 "astair.mods.gz",
                 "astair.stats")
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.fname.split('_')[1]],
        method       = lambda wildcards: config["astair"]["method"][wildcards.fname.split('_')[0]],
        extra_params = (config["astair"]["count"]["extra_params"]
                        if config["astair"]["count"]["extra_params"] else ""),
        pattern      = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        astair call \
            -i {wildcards.fname}{params.pattern}.bam \
            --reference {params.ref} \
            -co CpG -ni . -d astair -t {threads} \
            --method {params.method} {params.extra_params}
        cd astair
        mv \
            {wildcards.fname}{params.pattern}_{params.method}_CpG.mods \
            {wildcards.fname}{params.pattern}.astair.mods
        pigz --best {wildcards.fname}{params.pattern}.astair.mods
        mv \
            {wildcards.fname}{params.pattern}_{params.method}_CpG.stats \
            {wildcards.fname}{params.pattern}.astair.stats
        """
