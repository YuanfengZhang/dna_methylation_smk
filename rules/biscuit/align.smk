configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule biscuit_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz"
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/biscuit/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/biscuit/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/biscuit/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["biscuit"][wildcards.BaseName.split('_')[1]],
        extra_params = config["biscuit"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        cd result/{wildcards.BaseName}/{wildcards.AlignParentDir}
        mkdir -p biscuit
        biscuit align \\
            -@ {threads} \\
            -C -V -Y {params.extra_params} \\
            {params.ref} \\
            {wildcards.BaseName}.R1.fq.gz \\
            {wildcards.BaseName}.R2.fq.gz |\\
        mbuffer -m 4G -q |\\
        samtools sort - |\\
        samtools addreplacerg - \\
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            --output-fmt bam,level=9 \\
            -@ {threads} -o biscuit/{wildcards.BaseName}.bam

        samtools index -@ {threads} biscuit/{wildcards.BaseName}.bam || echo "suppress non-zero exit"
        """)
