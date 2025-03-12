configfile: "config/runtime_config.yaml"

rule biscuit_count:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        ("result/{fname}/{trimmer}/{aligner}/{deduper}/"
         "biscuit/{fname}{bqsr}epibed.gz")
    params:
        ref          = lambda wildcards: config["ref"]["biscuit"][wildcards.fname.split('_')[1]],
        pileup_params = (config["biscuit"]["pileup"]["extra_params"]
                        if config["biscuit"]["pileup"]["extra_params"] else ""),
        vcf2bed_params = (config["biscuit"]["vcf2bed"]["extra_params"]
                        if config["biscuit"]["vcf2bed"]["extra_params"] else ""),
        pattern      = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        mkdir -p biscuit
        biscuit pileup \
            -@ {threads} -o biscuit/{wildcards.fname}{params.pattern}.vcf \
            {params.ref} {wildcards.fname}{params.pattern}.bam {params.pileup_params} |\
        bgzip -@ {threads} -o biscuit/{wildcards.fname}{params.pattern}.vcf.gz
        tabix -p vcf biscuit/{wildcards.fname}{params.pattern}.vcf.gz
        biscuit vcf2bed \
            -e {params.vcf2bed_params} biscuit/{wildcards.fname}{params.pattern}.vcf |\
        pigz --best > biscuit/{wildcards.fname}{params.pattern}.epibed.gz

        rm biscuit/{wildcards.fname}{params.pattern}.vcf.gz
        rm biscuit/{wildcards.fname}{params.pattern}.vcf.gz.tbi
        """

