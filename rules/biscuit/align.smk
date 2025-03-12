configfile: "config/runtime_config.yaml"

rule biscuit_align:
    input:
        r1           = "result/{fname}/{trimmer}/{fname}.R1.fq.gz",
        r2           = "result/{fname}/{trimmer}/{fname}.R2.fq.gz"
    output:
        bam          = "result/{fname}/{trimmer}/biscuit/{fname}.bam",
        bai          = "result/{fname}/{trimmer}/biscuit/{fname}.bam.bai"
    params:
        ref          = lambda wildcards: config["ref"]["biscuit"][wildcards.fname.split('_')[1]],
        extra_params = (config["biscuit"]["align"]["extra_params"]
                        if config["biscuit"]["align"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        LIB=$(echo "{wildcards.fname}" | cut -d'_' -f1)
        PLATFORM=$(echo "{wildcards.fname}" | cut -d'_' -f4)
        SAMPLE=$(echo "{wildcards.fname}" | cut -d'_' -f2-3)

        cd result/{wildcards.fname}/{wildcards.trimmer}
        mkdir -p biscuit
        biscuit align \
            -@ {threads} \
            -C -V -Y {params.extra_params} \
            {params.ref} {wildcards.fname}.R1.fq.gz {wildcards.fname}.R2.fq.gz |\
        mbuffer -m 8G |\
        samtools sort - |\
        samtools addreplacerg - \
            -r "@RG\\tID:{wildcards.fname}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}\" \
            --output-fmt bam,level=9 \
            -@ {threads} -o biscuit/{wildcards.fname}.bam

        samtools index -@ {threads} biscuit/{wildcards.fname}.bam || echo "suppress non-zero exit"
        """
