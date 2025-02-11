configfile: "config/runtime_config.yaml"

rule qualimap_bqsr:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}.bqsr.bam"
    output:
        multiext("result/{fname}/{trimmer}/{aligner}/"
                 "{deduper}/qualimap_bqsr/",
                 "{fname}.bqsr.qualimap.pdf",
                 "qualimapReport.html")
    params:
        extra_params = (config["qualimap"]["extra_params"]
                        if config["qualimap"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        qualimap bamqc \
            -bam {wildcards.fname}.bam \
            -outdir . \
            -outfile {wildcards.fname}.qualimap.pdf \
            -outformat PDF:HTML \
            -nt {threads} {params.extra_params}
        if [ -d "qualimap_bqsr" ]; then
            rm -rf qualimap_bqsr
        fi
        mv {wildcards.fname}{params.pattern}.bqsr_stats qualimap_bqsr
        """
