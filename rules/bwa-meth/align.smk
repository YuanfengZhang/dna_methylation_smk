configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule bwameth_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz"
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/bwa-meth/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/bwa-meth/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/bwa-meth/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-meth"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bwa-meth"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)
        bwameth.py \\
            --read-group "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            --threads {threads} \\
            --reference {params.ref} {params.extra_params} \\
            {input.r1} {input.r2} |\\
        mbuffer -m 4G -q |\\
        samtools sort \\
            -O bam,level=9 -@ {threads} \\
            -o {output.bam} -

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
