configfile: "config/runtime_config.yaml"
from textwrap import dedent

# We have to add fake SAM @RG headers to the BWA-MEME output
# so that the GATK MarkDuplicates can run

rule bwa_meme:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/bwa-meme/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/bwa-meme/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/bwa-meme/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-meme"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bwa-meme"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)
        bwa-meme mem \\
            -7 -Y -K 100000000 -t {threads} \\
            -R "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            {params.ref} {input.r1} {input.r2} |\\
        mbuffer -m 4G -q |\\
        samtools sort \\
            --output-fmt bam,level=9 \\
            -@ {threads} - >{output.bam}
        
        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
