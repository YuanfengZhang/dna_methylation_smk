configfile: "config/runtime_config.yaml"

rule qualimap:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}.bam"
    output:
        multiext("result/{fname}/{trimmer}/{aligner}/"
                 "{deduper}/qualimap/",
                 "{fname}.qualimap.pdf",
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
        if [ -d "qualimap" ]; then
            rm -rf qualimap
        fi
        mv {wildcards.fname}_stats qualimap
        """
