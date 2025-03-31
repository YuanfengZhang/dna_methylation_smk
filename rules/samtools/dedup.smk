configfile: "config/runtime_config.yaml"

rule samtools_dedup:
    input:
        "result/{BaseName}/{DedupParentDir}/{BaseName}.bam"
    output:
        bam            = "result/{BaseName}/{DedupParentDir}/samtools/{BaseName}.bam",
        bai            = "result/{BaseName}/{DedupParentDir}/samtools/{BaseName}.bam.bai",
        stats          = "result/{BaseName}/{DedupParentDir}/samtools/{BaseName}.dup.stats"
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
