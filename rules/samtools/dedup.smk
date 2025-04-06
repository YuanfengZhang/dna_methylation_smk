configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule samtools_dedup:
    input:
        "result/{BaseName}/{DedupParentDir}/{BaseName}.bam"
    output:
        bam            = "result/{BaseName}/{DedupParentDir}/samtools/{BaseName}.bam",
        bai            = "result/{BaseName}/{DedupParentDir}/samtools/{BaseName}.bam.bai",
        stats          = "result/{BaseName}/{DedupParentDir}/samtools/{BaseName}.dup.stats"
    params:
        fixmate_params = config["samtools"]["fixmate"]["extra_params"] or "",
        mkdup_params   = config["samtools"]["mkdup"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        samtools sort \
            -@ {threads} -n \
            -m 4G {input} |\
        samtools fixmate \
            -cmM -@ {threads} {params.fixmate_params} \
            - - |\
        samtools sort \
            -@ {threads} -m 4G |\
        samtools markdup \
            -@ {threads} -f {output.stats} {params.mkdup_params} \
            - {output.bam}
        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
