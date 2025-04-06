configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule biscuit_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        "result/{BaseName}/{CountParentDir}/biscuit/{BaseName}.epibed.gz"
    params:
        ref            = lambda wildcards: config["ref"]["biscuit"][wildcards.BaseName.split('_')[1]],
        pileup_params  = config["biscuit"]["pileup"]["extra_params"] or "",
        vcf2bed_params = config["biscuit"]["vcf2bed"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.CountParentDir}
        mkdir -p biscuit
        biscuit pileup \
            -@ {threads} -o biscuit/{wildcards.BaseName}.vcf \
            {params.ref} {wildcards.BaseName}.bam {params.pileup_params} |\
        bgzip -@ {threads} -o biscuit/{wildcards.BaseName}.vcf.gz
        tabix -p vcf biscuit/{wildcards.BaseName}.vcf.gz
        biscuit vcf2bed \
            -e {params.vcf2bed_params} biscuit/{wildcards.BaseName}.vcf |\
        pigz --best > biscuit/{wildcards.BaseName}.epibed.gz

        rm biscuit/{wildcards.BaseName}.vcf.gz
        rm biscuit/{wildcards.BaseName}.vcf.gz.tbi
        """)

