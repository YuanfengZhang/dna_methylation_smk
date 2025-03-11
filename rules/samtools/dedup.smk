configfile: "config/runtime_config.yaml"

rule samtools_dedup:
    input:
        "result/{fname}/{trimmer}/{aligner}/{fname}.bam"
    output:
        bam            = "result/{fname}/{trimmer}/{aligner}/samtools/{fname}.bam",
        bai            = "result/{fname}/{trimmer}/{aligner}/samtools/{fname}.bam.bai",
        stats          = "result/{fname}/{trimmer}/{aligner}/samtools/{fname}.dup.stats"
    params:
        fixmate_params = (config["samtools"]["fixmate"]["extra_params"]
                          if config["samtools"]["fixmate"]["extra_params"] else ""),
        mkdup_params   = (config["samtools"]["mkdup"]["extra_params"]
                          if config["samtools"]["mkdup"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        samtools sort \
            -@ {threads} -n \
            -m 1G {input} |\
        samtools fixmate \
            -cmM -@ {threads} {params.fixmate_params} \
            - - |\
        samtools sort \
            -@ {threads} -m 1G |\
        samtools markdup \
            -@ {threads} -f {output.stats} {params.mkdup_params} \
            - {output.bam}
        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """
