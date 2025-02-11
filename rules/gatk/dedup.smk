configfile: "config/runtime_config.yaml"

rule gatk_dedup:
    input:
        "result/{fname}/{trimmer}/{aligner}/{fname}.bam"
    output:
        bam          = "result/{fname}/{trimmer}/{aligner}/gatk/{fname}.bam",
        bai          = "result/{fname}/{trimmer}/{aligner}/gatk/{fname}.bam.bai",
        metrics      = "result/{fname}/{trimmer}/{aligner}/gatk/{fname}.dedup.metrics"
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
            --ASSUME_SORT_ORDER coordinate {params.extra_params}

        samtools index -@ {threads} {output.bam}
        """
