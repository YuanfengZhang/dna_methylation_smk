configfile: "config/runtime_config.yaml"
from textwrap import dedent

# https://github.com/bioinformatics-polito/BioSeqZip

rule BioSeqZip:
    input:
        r1           = "result/{BaseName}/{DedupParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{DedupParentDir}/{BaseName}.R2.fq.gz"
    output:
        r1           = "result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.R2.fq.gz",
        report       = "result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.bioseqzip.report"
    benchmark:
        "result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.dedup.benchmark"
    params:
        extra_params = lambda wildcards: config["bio-seq-zip"]["dedup"]["extra_params"] or ""
    threads: 8
    conda:
        "../biscuit/conda.yaml"
    shell:
        dedent("""
        export tmp_dir="result/{wildcards.BaseName}/{wildcards.DedupParentDir}/bio-seq-zip"
        mkdir -p ${{tmp_dir}}
        resources/BioSeqZip/build/apps/bioseqzip_collapse \\
            -i {input.r1} \\
            -p {input.r2} \\
            -f fastq \\
            -o ${{tmp_dir}} \\
            -b {wildcards.BaseName} \\
            --csv-report -v 4 {params.extra_params}
        mv \\
          ${{tmp_dir}}/{wildcards.BaseName}_1.fq \\
          ${{tmp_dir}}/{wildcards.BaseName}.R1.fq
        mv \\
          ${{tmp_dir}}/{wildcards.BaseName}_2.fq \\
          ${{tmp_dir}}/{wildcards.BaseName}.R2.fq
        mv \\
          ${{tmp_dir}}/{wildcards.BaseName}.csv \\
          ${{tmp_dir}}/{wildcards.BaseName}.bioseqzip.report
        pigz --best -p {threads} ${{tmp_dir}}/{wildcards.BaseName}.R1.fq
        pigz --best -p {threads} ${{tmp_dir}}/{wildcards.BaseName}.R2.fq
        """)
