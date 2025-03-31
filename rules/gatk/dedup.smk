configfile: "config/runtime_config.yaml"

rule gatk_dedup:
    input:
        "result/{BaseName}/{DedupParentDir}/{BaseName}.bam"
    output:
        bam          = "result/{BaseName}/{DedupParentDir}/gatk-dedup/{BaseName}.bam",
        bai          = "result/{BaseName}/{DedupParentDir}/gatk-dedup/{BaseName}.bam.bai",
        metrics      = "result/{BaseName}/{DedupParentDir}/gatk-dedup/{BaseName}.dedup.metrics"
    params:
        extra_params = (config["gatk"]["dedup"]["extra_params"]
                        if config["gatk"]["dedup"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        gatk \
            --java-options "-Xmx20g -XX:ParallelGCThreads={threads}" \
            MarkDuplicates \
            -I {input} \
            -O {output.bam} \
            -M {output.metrics} \
            --ASSUME_SORT_ORDER coordinate \
            {params.extra_params}

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """
