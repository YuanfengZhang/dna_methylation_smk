configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
bowtie2-build --large-index --threads 64 BL.fa BL
bowtie2 -x genome \
    -1 reads_1.fq \
    -2 read2_2.fq \
    -S output.sam \
    -q --sensitive \
    --rg-id ID:1 \
    --rg SM:sample1 \
    --rg PL:ILLUMINA \
    --rg PU:unit1 \
    --rg LB:lib1 \
    -t -p 16
"""


rule bowtie2_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/bowtie2/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/bowtie2/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/bowtie2/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bowtie2"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bowtie2"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d'_' -f1)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d'_' -f4)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        bowtie2 \
            -x {params.ref} -1 {input.r1} -2 {input.r2} \
            --rg-id ID:{wildcards.BaseName} \
            --rg SM:${{SAMPLE}} \
            --rg PL:${{PLATFORM}}\
            --rg LB:${{LIB}} \
            -p {threads} {params.extra_params} |\
        samtools sort \
            --output-fmt bam,level=9 \
            -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
