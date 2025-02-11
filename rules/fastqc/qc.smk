configfile: "config/runtime_config.yaml"

rule fastqc:
    input:
        r1           = "result/{fname}/{trimmer}/{fname}.R1.fq.gz",
        r2           = "result/{fname}/{trimmer}/{fname}.R2.fq.gz"
    output:
        r1_html      = "result/{fname}/{trimmer}/{fname}.R1.fastqc.html",
        r2_html      = "result/{fname}/{trimmer}/{fname}.R2.fastqc.html"
    params:
        extra_params = config["fastqc"]["extra_params"]
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        fastqc \
            {params.extra_params} \
            --threads {threads} \
            -o result/{wildcards.fname}/{wildcards.trimmer} \
            {input.r1} {input.r2}
        """
