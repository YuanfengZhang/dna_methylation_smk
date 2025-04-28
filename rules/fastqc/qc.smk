configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule fastqc:
    input:
        r1           = "result/{BaseName}/{ReportParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{ReportParentDir}/{BaseName}.R2.fq.gz"
    output:
        r1_html      = "result/{BaseName}/{ReportParentDir}/{BaseName}.R1.fastqc.html",
        r2_html      = "result/{BaseName}/{ReportParentDir}/{BaseName}.R2.fastqc.html"
    benchmark:
        "result/{BaseName}/{ReportParentDir}/{BaseName}.fastqc.benchmark"
    params:
        extra_params = config["fastqc"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        fastqc {params.extra_params} \\
            --threads {threads} \\
            -o result/{wildcards.BaseName}/{wildcards.ReportParentDir} \\
            {input.r1} {input.r2}
        """)