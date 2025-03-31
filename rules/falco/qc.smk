configfile: "config/runtime_config.yaml"

rule falco:
    input:
        r1            = "result/{BaseName}/{ReportParentDir}/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{ReportParentDir}/{BaseName}.R2.fq.gz"
    output:
        multiext("result/{BaseName}/{ReportParentDir}/{BaseName}",
                 ".R1.falco.html", ".R2.falco.html",
                 ".R1.falco.data.txt", ".R2.falco.data.txt",
                 ".R1.falco.summary", ".R2.falco.summary")
    params:
        extra_params  = (config["falco"]["extra_params"]
                         if config["falco"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        falco \
            --threads {threads} \
            {input.r1} {input.r2} \
            -bisulfite {params.extra_params}
        cd result/{wildcards.BaseName}/{wildcards.ReportParentDir}

        mv {wildcards.BaseName}.R1.fq.gz_fastqc_report.html  {wildcards.BaseName}.R1.falco.html
        mv {wildcards.BaseName}.R2.fq.gz_fastqc_report.html {wildcards.BaseName}.R2.falco.html
        mv {wildcards.BaseName}.R1.fq.gz_fastqc_data.txt {wildcards.BaseName}.R1.falco.data.txt
        mv {wildcards.BaseName}.R2.fq.gz_fastqc_data.txt {wildcards.BaseName}.R2.falco.data.txt
        mv {wildcards.BaseName}.R1.fq.gz_summary.txt {wildcards.BaseName}.R1.falco.summary
        mv {wildcards.BaseName}.R2.fq.gz_summary.txt {wildcards.BaseName}.R2.falco.summary
        """
