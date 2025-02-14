#! DO NOT USE THIS! USE BWA-MEME INSTEAD
configfile: "config/runtime_config.yaml"

rule astair_align:
    input:
        r1           = "result/{fname}/{trimmer}/{fname}.R1.fq.gz",
        r2           = "result/{fname}/{trimmer}/{fname}.R1.fq.gz"
    output:
        "result/{fname}/{trimmer}/astair/{fname}.bam"
    params:
        ref          = lambda wildcards: config["ref"]["astair"][wildcards.fname.split('_')[1]],
        method       = lambda wildcards: config["astair"]["method"][wildcards.fname.split('_')[0]],
        extra_params = (config["astair"]["align"]["extra_params"]
                        if config["astair"]["align"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        astair align \
            --fq1 {input.r1} --fq2 {input.r2} \
            --reference {params.ref} \
            -d result/{wildcards.fname}/{wildcards.trimmer}/astair \
            -n {wildcards.fname} \
            -t {threads} --method {params.method} {params.extra_params}
        cd result/{wildcards.fname}/{wildcards.trimmer}/astair
        samtools sort -O bam,level=9 \
            -o {wildcards.fname}.bam \
            {wildcards.fname}_{params.method}.bam
        rm {wildcards.fname}_{params.method}.bam
        samtools index -@ {threads} {wildcards.fname}.bam || echo "suppress non-zero exit"
        """
