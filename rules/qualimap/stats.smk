configfile: "config/runtime_config.yaml"

rule qualimap:
    input:
        "result/{BaseName}/{QualimapParentDir}/{BaseName}.bam"
    output:
        multiext("result/{BaseName}/{QualimapParentDir}/qualimap/",
                 "{BaseName}.qualimap.pdf",
                 "qualimapReport.html")
    params:
        extra_params = (config["qualimap"]["extra_params"]
                        if config["qualimap"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.BaseName}/{wildcards.QualimapParentDir}
        if [ -d "qualimap" ]; then
            rm -rf qualimap
        fi

        qualimap bamqc \
            -bam {wildcards.BaseName}.bam \
            -outdir qualimap \
            -outfile {wildcards.BaseName}.qualimap.pdf \
            -outformat PDF:HTML \
            -nt {threads} {params.extra_params}
        """
