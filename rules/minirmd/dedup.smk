configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule minirmd_dedup:
    input:
        r1           = "result/{BaseName}/{DedupParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{DedupParentDir}/{BaseName}.R2.fq.gz"
    output:
        tmp1         = temp("result/{BaseName}/{DedupParentDir}/minirmd/{BaseName}_1"),
        tmp2         = temp("result/{BaseName}/{DedupParentDir}/minirmd/{BaseName}_2"),
        r1           = "result/{BaseName}/{DedupParentDir}/minirmd/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{DedupParentDir}/minirmd/{BaseName}.R2.fq.gz"
    benchmark:
        "result/{BaseName}/{DedupParentDir}/minirmd/{BaseName}.dedup.benchmark"
    params:
        mismatch_num = lambda wildcards: config["minirmd"]["dedup"]["mismatch_num"],
        extra_params = lambda wildcards: config["minirmd"]["dedup"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml" 
    shell:
        dedent("""
        minirmd \\
            -i {input.r1} \\
            -f {input.r2} \\
            -o "result/{wildcards.BaseName}/{wildcards.DedupParentDir}/minirmd/{wildcards.BaseName}" \\
            -d {params.mismatch_num} \\
            -t {threads} {params.extra_params}
        gzip -c {output.tmp1} > {output.r1}
        gzip -c {output.tmp2} > {output.r2}
        """)
