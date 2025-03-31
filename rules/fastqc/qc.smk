configfile: "config/runtime_config.yaml"

rule fastqc:
    input:
        r1           = "result/{BaseName}/{ReportParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{ReportParentDir}/{BaseName}.R2.fq.gz"
    output:
        r1_html      = "result/{BaseName}/{ReportParentDir}/{BaseName}.R1.fastqc.html",
        r2_html      = "result/{BaseName}/{ReportParentDir}/{BaseName}.R2.fastqc.html"
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
            -o result/{wildcards.BaseName}/{wildcards.ReportParentDir} \
            {input.r1} {input.r2}
        """
