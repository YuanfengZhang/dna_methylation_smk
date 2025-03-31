configfile: "config/runtime_config.yaml"

rule bwameth_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz"
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/bwa-meth/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/bwa-meth/{BaseName}.bam.bai"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-meth"][wildcards.BaseName.split('_')[1]],
        extra_params = (config["bwa-meth"]["extra_params"]
                        if config["bwa-meth"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        LIB=$(echo "{wildcards.BaseName}" | cut -d'_' -f1)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d'_' -f4)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)
        bwameth.py \
            --read-group "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}\" \
            --threads {threads} \
            --reference {params.ref} {params.extra_params} \
            {input.r1} {input.r2} |\
        mbuffer -m 1G |\
        samtools sort \
            -O bam,level=9 -@ {threads} \
            -o {output.bam} -

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """
