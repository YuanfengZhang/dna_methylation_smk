configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
https://github.com/bioinformatics-polito/BioSeqZip

multi CPUs do not help:
#######################################################
32 cores + 12G

s       h:m:s   max_rss max_vms max_uss max_pss 
257.8878        0:04:17 3565.51 5633.54 3572.53

io_in   io_out  mean_load cpu_time
3572.94 7320.85 11264.01 1074.83 3416.83
#######################################################
1 core + 4G

s       h:m:s   max_rss max_vms max_uss max_pss
291.9358        0:04:51 3390.46 3402.65 3384.05 3384.46

io_in   io_out  mean_load       cpu_time
8291.38 11212.11        924.13  3175.97
"""

rule BioSeqZip:
    input:
        r1           = "result/{BaseName}/{DedupParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{DedupParentDir}/{BaseName}.R2.fq.gz"
    output:
        r1           = temp("result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.R1.fq"),
        r2           = temp("result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.R2.fq"),
        report       = "result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.bioseqzip.report"
    benchmark:
        "result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.dedup.benchmark"
    params:
        extra_params = config["bio-seq-zip"]["dedup"]["extra_params"] or ""
    threads: 1  # threads won't be fully utilized
    resources:
        mem_mb = 8192
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
            -t {threads} \\
            --tmp-directory ${{tmp_dir}} \\
            --output-directory ${{tmp_dir}} \\
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
        """)

rule pigz_compress:
    input:
        r1 = "result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.R1.fq",
        r2 = "result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.R2.fq"
    output:
        r1 = "result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.R1.fq.gz",
        r2 = "result/{BaseName}/{DedupParentDir}/bio-seq-zip/{BaseName}.R2.fq.gz"
    threads: 8
    conda:
        "../biscuit/conda.yaml"
    shell:
        dedent("""
        pigz --best -p {threads} -k {input.r1}
        pigz --best -p {threads} -k {input.r2}
        """)
