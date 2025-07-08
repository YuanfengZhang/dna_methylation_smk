configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule ngs_reads_treatment:
    input:
        r1            = "result/{BaseName}/{DedupParentDir}/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{DedupParentDir}/{BaseName}.R2.fq.gz"
    output:
        tmp_inr1      = temp("result/{BaseName}/{DedupParentDir}/ngs-reads-treatment/{BaseName}.R1.fq"),
        tmp_inr2      = temp("result/{BaseName}/{DedupParentDir}/ngs-reads-treatment/{BaseName}.R2.fq"),
        r1            = "result/{BaseName}/{DedupParentDir}/ngs-reads-treatment/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{DedupParentDir}/ngs-reads-treatment/{BaseName}.R2.fq.gz"
    benchmark:
        "result/{BaseName}/{DedupParentDir}/ngs-reads-treatment/{BaseName}.dedup.benchmark"
    params:
        max_mem       = config["ngs-reads-treatment"]["dedup"]["max_mem"],
        extra_params  = lambda wildcards: config["ngs-reads-treatment"]["dedup"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        export tmp_dir="result/{wildcards.BaseName}/{wildcards.DedupParentDir}/ngs-reads-treatment"
        mkdir -p ${{tmp_dir}}
        pigz -dc {input.r1} > {output.tmp_inr1}
        pigz -dc {input.r2} > {output.tmp_inr2}
        java \\
          -Xmx{params.max_mem}G \\
          -jar resources/ngsreadstreatment/NgsReadsTreatment_v1.3.jar \\
          {output.tmp_inr1} {output.tmp_inr2} \\
          {threads} {params.extra_params}

        mv \\
          ${{tmp_dir}}/{wildcards.BaseName}_1_trated.fastq \\
          ${{tmp_dir}}/{wildcards.BaseName}.R1.fq
        
        mv \\
          ${{tmp_dir}}/{wildcards.BaseName}_2_trated.fastq \\
          ${{tmp_dir}}/{wildcards.BaseName}.R2.fq

        pigz --best -p {threads} ${{tmp_dir}}/{wildcards.BaseName}.R1.fq
        pigz --best -p {threads} ${{tmp_dir}}/{wildcards.BaseName}.R2.fq
        """)
