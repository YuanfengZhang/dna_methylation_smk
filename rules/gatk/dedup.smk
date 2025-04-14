configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule gatk_dedup:
    input:
        "result/{BaseName}/{DedupParentDir}/{BaseName}.bam"
    output:
        bam          = "result/{BaseName}/{DedupParentDir}/gatk-dedup/{BaseName}.bam",
        bai          = "result/{BaseName}/{DedupParentDir}/gatk-dedup/{BaseName}.bam.bai",
        metrics      = "result/{BaseName}/{DedupParentDir}/gatk-dedup/{BaseName}.dedup.metrics"
    benchmark:
        "result/{BaseName}/{DedupParentDir}/gatk-dedup/{BaseName}.dedup.benchmark"
    params:
        extra_params = config["gatk"]["dedup"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        gatk \
            --java-options "-Xmx20g -XX:ParallelGCThreads={threads}" \
            MarkDuplicates \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} \
            --ASSUME_SORT_ORDER coordinate \
            {params.extra_params}

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
