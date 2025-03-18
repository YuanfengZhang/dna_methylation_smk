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
        extra_params = (config["bsgenova"]["count"]["extra_params"]
                        if config["bsgenova"]["count"]["extra_params"] else ""),
        pattern      = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        bsgenova_dir="result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}"
        python resources/bsgenova/bsextractor.py \
            -b {wildcards.fname}{params.pattern}.bam \
            -g {params.ref} \
            --output-atcgmap ${{bsgenova_dir}}/{wildcards.fname}{params.pattern}.bsgenova.ATCGmap.gz \
            --output-cgmap ${{bsgenova_dir}}/{wildcards.fname}{params.pattern}.bsgenova.CGmap.gz \
            --output-bed ${{bsgenova_dir}}/{wildcards.fname}{params.pattern}.bsgenova.bed.gz \
            --threads {threads} {params.extra_params} 
        """
