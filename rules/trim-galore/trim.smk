configfile: "config/runtime_config.yaml"

rule trim_galore:
    input:
        r1           = (f"{config['input_dir']}"
                        "/{fname}.R1.fq.gz"),
        r2           = (f"{config['input_dir']}"
                        "/{fname}.R2.fq.gz")
    output:
        r1           = "result/{fname}/trim-galore/{fname}.R1.fq.gz",
        r2           = "result/{fname}/trim-galore/{fname}.R2.fq.gz",
        r1_report    = "result/{fname}/trim-galore/{fname}.R1.report",
        r2_report    = "result/{fname}/trim-galore/{fname}.R2.report"
    params:
        extra_params = config["trim-galore"]["extra_params"]
    threads: 8
    conda:
        "rule/trim-galore/conda.yaml"
    shell:
        """
        trim_galore \
            {params.extra_params} \
            -j {threads} \
            --output_dir result/{wildcards.fname}/trim-galore \
            {input.r1} {input.r2}
        mv \
            result/{wildcards.fname}/trim-galore/{wildcards.fname}.R1_trimmed.fq.gz \
            {output.r1}
        mv \
            result/{wildcards.fname}/trim-galore/{wildcards.fname}.R2_trimmed.fq.gz \
            {output.r2}
        mv \
            result/{wildcards.fname}/trim-galore/{wildcards.fname}.R1.fq.gz_trimming_report.txt \
            {output.r1_report}
        mv \
            result/{wildcards.fname}/trim-galore/{wildcards.fname}.R2.fq.gz_trimming_report.txt \
            {output.r2_report}
        """