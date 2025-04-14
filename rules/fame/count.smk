configfile: "config/runtime_config.yaml"
from textwrap import dedent


# to build the index:
# resources/FAME/FAME --genome ref.fa --store_index ref.fa.fame

rule fq2tsv:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz"
    output:
        'result/{BaseName}/{AlignParentDir}/fame/{BaseName}.bedgraph.zst'
    benchmark:
        "result/{BaseName}/{AlignParentDir}/fame/{BaseName}.count.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["fame"][wildcards.BaseName.split('_')[1]],
        extra_params = config["fame"]["extra_params"] or ""
    threads: 8
    shell:
        dedent("""
        resources/FAME/FAME \
            -r1 {input.r1} -r2 {input.r2} \
            --load_index {params.ref} {params.extra_params} \
            -o result/{wildcards.BaseName}/{wildcards.AlignParentDir}/fame/{wildcards.BaseName}
        POLARS_MAX_THREADS={threads} python rules/fame/format.py \
            --input result/{wildcards.BaseName}/{wildcards.AlignParentDir}/fame/{wildcards.BaseName}_cpg.tsv \
            --output {output} \
            --threads {threads}
        #rm result/{wildcards.BaseName}/{wildcards.AlignParentDir}/fame/{wildcards.BaseName}_cpg.tsv
        """)
