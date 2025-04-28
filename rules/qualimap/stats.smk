configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule qualimap:
    input:
        "result/{BaseName}/{QualimapParentDir}/{BaseName}.bam"
    output:
        multiext("result/{BaseName}/{QualimapParentDir}/qualimap/",
                 "{BaseName}.qualimap.pdf",
                 "qualimapReport.html")
    benchmark:
        "result/{BaseName}/{QualimapParentDir}/qualimap/{BaseName}.qualimap.benchmark"
    params:
        extra_params = config["qualimap"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.QualimapParentDir}
        if [ -d "qualimap" ]; then
            rm -rf qualimap
        fi

        qualimap bamqc \\
            -bam {wildcards.BaseName}.bam \\
            -outdir qualimap \\
            -outfile {wildcards.BaseName}.qualimap.pdf \\
            -outformat PDF:HTML \\
            -nt {threads} {params.extra_params}
        """)
