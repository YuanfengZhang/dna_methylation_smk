configfile: "config/runtime_config.yaml"

rule fastp:
    input:
        r1           = (f"{config.input_dir}"
                        "/{fname}.R1.fq.gz"),
        r2           = (f"{config.input_dir}"
                        "/{fname}.R2.fq.gz")
    output:
        r1           = "result/{fname}/fastp/{fname}.R1.fq.gz",
        r2           = "result/{fname}/fastp/{fname}.R2.fq.gz",
        qc_html      = "result/{fname}/fastp/{fname}.fastp.html",
        qc_json      = "result/{fname}/fastp/{fname}.fastp.json"
    params:
        extra_params = (config["fastp"]["extra_params"]
                        if config["fastp"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        fastp \
            -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            -h {output.qc_html} -j {output.qc_json} \
            -w {threads} \
            {params.extra_params}
        """
