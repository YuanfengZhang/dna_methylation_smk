configfile: "config/runtime_config.yaml"

rule all:
    input:
        "result/PS_HF_1_MX/trim-galore/PS_HF_1_MX.R1.fq.gz",
        "result/PS_HF_1_MX/trim-galore/PS_HF_1_MX.R2.fq.gz"

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
    log:
        stdout       = "result/{fname}/trim-galore/{fname}.stdout",
        stderr       = "result/{fname}/trim-galore/{fname}.stderr"
    params:
        extra_params = (config["trim-galore"]["extra_params"]
                        if config["trim-galore"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        trim_galore \
            -j {threads} --paired \
            --output_dir result/{wildcards.fname}/trim-galore \
            {input.r1} {input.r2} {params.extra_params} 2>{log.stderr} 1>{log.stdout}
        mv \
            result/{wildcards.fname}/trim-galore/{wildcards.fname}.R1_val_1.fq.gz \
            {output.r1}
        mv \
            result/{wildcards.fname}/trim-galore/{wildcards.fname}.R2_val_2.fq.gz \
            {output.r2}
        mv \
            result/{wildcards.fname}/trim-galore/{wildcards.fname}.R1.fq.gz_trimming_report.txt \
            {output.r1_report}
        mv \
            result/{wildcards.fname}/trim-galore/{wildcards.fname}.R2.fq.gz_trimming_report.txt \
            {output.r2_report}
        """
