configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule bsgenova_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        atcgmap      = "result/{BaseName}/{CountParentDir}/bsgenova/{BaseName}.ATCGmap.gz",
        cgmap        = "result/{BaseName}/{CountParentDir}/bsgenova/{BaseName}.CGmap.gz",
        bedgz        = "result/{BaseName}/{CountParentDir}/bsgenova/{BaseName}.bed.gz"
    benchmark:
        "result/{BaseName}/{CountParentDir}/bsgenova/{BaseName}.count.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bsgenova"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bsgenova"]["count"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        python resources/bsgenova/bsextractor.py \\
            -b {input} \\
            -g {params.ref} \\
            --output-atcgmap {output.atcgmap} \\
            --output-cgmap {output.cgmap} \\
            --output-bed {output.bedgz} \\
            --threads {threads} {params.extra_params} 
        """)
