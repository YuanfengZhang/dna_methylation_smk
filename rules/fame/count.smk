configfile: "config/runtime_config.yaml"

rule fq2tsv:
    input:
        r1           = "result/{fname}/{trimmer}/{fname}.R1.fq.gz",
        r2           = "result/{fname}/{trimmer}/{fname}.R2.fq.gz"
    output:
        'result/{fname}/{trimmer}/fame/{fname}.fame.bedgraph.gz'
    params:
        ref          = lambda wildcards: config["fame"]["ref"][wildcards.fname.split('_')[1]],
        extra_params = (config["fame"]["extra_params"]
                        if config["fame"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        resources/FAME/FAME \
            -r1 {input.r1} -r2 {input.r2} \
            --load_index {params.ref} {params.extra_param} \
            -o result/{wildcards.fname}/fame/{wildcards.fname}
        python rules/fame/format.py \
            --input result/{wildcards.fname}/{wildcards.trimmer}/fame/{wildcards.fname}_cpg.tsv \
            --output {output} \
            --threads {threads}
        rm result/{wildcards.fname}/{wildcards.trimmer}/fame/{wildcards.fname}_cpg.tsv
        """
