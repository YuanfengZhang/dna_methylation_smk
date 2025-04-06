configfile: "config/runtime_config.yaml"
from textwrap import dedent


#! no --rg_tag in bismark, but samtools addreplacerg after it,
# so we can make sure the fake @RG is identical to the output from bwa series.
#! samtools sort | samtools addreplacerg is faster than samtools addreplacerg | samtools sort.

rule bismark_bowtie2:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz"
    output:
        multiext("result/{BaseName}/{AlignParentDir}/bismark-bowtie2/{BaseName}",
                 ".bam",
                 ".align.report")
    params:
       # The BaseName example: BS_D5_1_Biochain.
       # We split the BaseName by '_' and get the second element.
       # ! Watch out that the dir containing "Bisulfite_Genome"
       # ! should be provided instead of the "Bisulfite_Genome" itself.
        ref          = lambda wildcards: config["ref"]["bismark-bowtie2"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bismark"]["align_bowtie2"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        bismark \
            -p {threads} --bam \
            --output_dir result/{wildcards.BaseName}/{wildcards.AlignParentDir}/bismark-bowtie2 \
            --genome_folder {params.ref} \
            -1 {input.r1} -2 {input.r2} {params.extra_params}

        cd result/{wildcards.BaseName}/{wildcards.AlignParentDir}/bismark-bowtie2

        LIB=$(echo "{wildcards.BaseName}" | cut -d'_' -f1)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d'_' -f4)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        samtools sort \
            {wildcards.BaseName}.R1_bismark_bt2_pe.bam \
            -@ {threads} \
            -m 4G |\
        samtools addreplacerg - \
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \
            --output-fmt bam,level=9 \
            -@ {threads} -o {wildcards.BaseName}.bam

        samtools index -@ {threads} {wildcards.BaseName}.bam || echo "suppress non-zero exit"

        mv \
            {wildcards.BaseName}.R1_bismark_bt2_PE_report.txt \
            {wildcards.BaseName}.align.report

        rm {wildcards.BaseName}.R1_bismark_bt2_pe.bam
        """)

rule bismark_hisat2:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz"
    output:
        multiext("result/{BaseName}/{AlignParentDir}/bismark-hisat2/{BaseName}",
                 ".bam",
                 ".align.report")
    params:
       # The BaseName example: BS_D5_1_Biochain.
       # We split the BaseName by '_' and get the second element.
        ref          = lambda wildcards: config["ref"]["bismark-hisat2"][wildcards.BaseName.split('_')[1]],
        extra_params = config["bismark"]["align_hisat2"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        bismark \
            -p {threads} --bam --hisat2 \
            --output_dir result/{wildcards.BaseName}/{wildcards.AlignParentDir}/bismark-hisat2 \
            --genome_folder {params.ref} \
            -1 {input.r1} -2 {input.r2} {params.extra_params}

        cd result/{wildcards.BaseName}/{wildcards.AlignParentDir}/bismark-hisat2

        LIB=$(echo "{wildcards.BaseName}" | cut -d'_' -f1)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d'_' -f4)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        samtools sort \
            {wildcards.BaseName}.R1_bismark_bt2_pe.bam \
            -m 4G |\
        samtools addreplacerg - \
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \
            --output-fmt bam,level=9 \
            -@ {threads} -o {wildcards.BaseName}.bam

        samtools index -@ {threads} {wildcards.BaseName}.bam || echo "suppress non-zero exit"

        mv \
            {wildcards.BaseName}.R1_bismark_bt2_PE_report.txt \
            {wildcards.BaseName}.align.report

        rm {wildcards.BaseName}.R1_bismark_bt2_pe.bam
        """)
