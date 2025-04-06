configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule trim_galore:
    input:
        r1           = (f"{config['input_dir']}/"
                         "{BaseName}.R1.fq.gz"),
        r2           = (f"{config['input_dir']}/"
                         "{BaseName}.R2.fq.gz")
    output:
        r1           = "result/{BaseName}/trim-galore/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/trim-galore/{BaseName}.R2.fq.gz",
        r1_report    = "result/{BaseName}/trim-galore/{BaseName}.R1.report",
        r2_report    = "result/{BaseName}/trim-galore/{BaseName}.R2.report"
    params:
        extra_params = config["trim-galore"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        trim_galore \
            -j {threads} --paired \
            --output_dir result/{wildcards.BaseName}/trim-galore \
            {input.r1} {input.r2} {params.extra_params}
        mv \
            result/{wildcards.BaseName}/trim-galore/{wildcards.BaseName}.R1_val_1.fq.gz \
            {output.r1}
        mv \
            result/{wildcards.BaseName}/trim-galore/{wildcards.BaseName}.R2_val_2.fq.gz \
            {output.r2}
        mv \
            result/{wildcards.BaseName}/trim-galore/{wildcards.BaseName}.R1.fq.gz_trimming_report.txt \
            {output.r1_report}
        mv \
            result/{wildcards.BaseName}/trim-galore/{wildcards.BaseName}.R2.fq.gz_trimming_report.txt \
            {output.r2_report}
        """)
