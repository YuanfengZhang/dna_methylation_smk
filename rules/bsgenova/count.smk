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
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bsgenova"]["count"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        tmp_dir="result/{wildcards.BaseName}/{wildcards.CountParentDir}/bsgenova/tmp"
        mkdir -p ${{tmp_dir}}
        python resources/bsgenova/bsextractor.py \\
            -b {input} \\
            -g {params.ref} \\
            --output-atcgmap {output.atcgmap} \\
            --output-cgmap {output.cgmap} \\
            --output-bed {output.bedgz} \\
            --tmp-dir ${{tmp_dir}} \\
            --threads {threads} {params.extra_params} 
        rm -rf ${{tmp_dir}}
        """)
