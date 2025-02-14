configfile: "config/runtime_config.yaml"

rule bwameth_align:
    input:
        r1           = "result/{fname}/{trimmer}/{fname}.R1.fq.gz",
        r2           = "result/{fname}/{trimmer}/{fname}.R2.fq.gz"
    output:
        bam          = "result/{fname}/{trimmer}/bwa-meth/{fname}.bam",
        bai          = "result/{fname}/{trimmer}/bwa-meth/{fname}.bam.bai"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-meth"][wildcards.fname.split('_')[1]],
        extra_params = (config["bwa-meth"]["extra_params"]
                        if config["bwa-meth"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        bwameth.py \
            --threads {threads} \
            --reference {params.ref} {params.extra_params} \
            {input.r1} {input.r2} |\
        mbuffer -m 1G |\
        samtools sort -m 1G \
            -O bam,level=9 -@ {threads} \
            -o {output.bam} -

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """
