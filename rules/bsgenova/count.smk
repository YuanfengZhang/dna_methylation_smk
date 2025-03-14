configfile: "config/runtime_config.yaml"

rule bsgenova_count:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        multiext("result/{fname}/{trimmer}/{aligner}/{deduper}"
                 "/bsgenova/{fname}{bqsr}",
                 "bsgenova.ATCGmap.gz",
                 "bsgenova.CGmap.gz",
                 "bsgenova.bed.gz")
    params:
        ref          = lambda wildcards: config["ref"]["bsgenova"][wildcards.fname.split('_')[1]],
        count_params = (config["bsgenova"]["count"]["extra_params"]
                        if config["bsgenova"]["count"]["extra_params"] else ""),
        pattern      = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        python ../../../../../resources/bsgenova/bsextractor.py \
            -b {wildcards.fname}{params.pattern}.bam \
            -g {params.ref} \
            --output-atcgmap bsgenova/{wildcards.fname}{params.pattern}.bsgenova.ATCGmap.gz \
            --output-cgmap bsgenova/{wildcards.fname}{params.pattern}.bsgenova.CGmap.gz \
            --output-bed bsgenova/{wildcards.fname}{params.pattern}.bsgenova.bed.gz \
            --threads {threads} {params.count_params} 
        """

