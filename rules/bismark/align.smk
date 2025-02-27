configfile: "config/runtime_config.yaml"

#! no --rg_tag in bismark, but samtools addreplacerg after it,
# so we can make sure the fake @RG is identical to the output from bwa series.
#! samtools sort | samtools addreplacerg is faster than samtools addreplacerg | samtools sort.

rule bismark_bowtie2:
    input:
        r1           = "result/{fname}/{trimmer}/{fname}.R1.fq.gz",
        r2           = "result/{fname}/{trimmer}/{fname}.R2.fq.gz"
    output:
        multiext("result/{fname}/{trimmer}/bismark-bowtie2/{fname}",
                 ".bam",
                 ".align.report")
    params:
       # The fname example: BS_D5_1_Biochain.
       # We split the fname by '_' and get the second element.
       # ! Watch out that the dir containing "Bisulfite_Genome"
       # ! should be provided instead of the "Bisulfite_Genome" itself.
        ref          = lambda wildcards: config["ref"]["bismark-bowtie2"][wildcards.fname.split('_')[1]],
        extra_params = (config["bismark"]["align_bowtie2"]["extra_params"]
                        if config["bismark"]["align_bowtie2"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        bismark \
            -p {threads} --bam \
            --output_dir result/{wildcards.fname}/{wildcards.trimmer}/bismark-bowtie2 \
            --genome_folder {params.ref} \
            -1 {input.r1} -2 {input.r2} {params.extra_params}

        cd result/{wildcards.fname}/{wildcards.trimmer}/bismark-bowtie2

        LIB=$(echo "{wildcards.fname}" | cut -d'_' -f1)
        PLATFORM=$(echo "{wildcards.fname}" | cut -d'_' -f4)
        SAMPLE=$(echo "{wildcards.fname}" | cut -d'_' -f2-3)

        samtools sort \
            {wildcards.fname}.R1_bismark_bt2_pe.bam \
            -m 1G |\
        samtools addreplacerg - \
            -r "@RG\\tID:{wildcards.fname}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}\" \
            --output-fmt bam,level=9 \
            -@ {threads} -o {wildcards.fname}.bam

        samtools index -@ {threads} {wildcards.fname}.bam || echo "suppress non-zero exit"

        mv \
            {wildcards.fname}.R1_bismark_bt2_PE_report.txt \
            {wildcards.fname}.align.report

        rm {wildcards.fname}.R1_bismark_bt2_pe.bam
        """

rule bismark_hisat2:
    input:
        r1           = "result/{fname}/{trimmer}/{fname}.R1.fq.gz",
        r2           = "result/{fname}/{trimmer}/{fname}.R2.fq.gz"
    output:
        multiext("result/{fname}/{trimmer}/bismark-hisat2/{fname}",
                 ".bam",
                 ".align.report")
    params:
       # The fname example: BS_D5_1_Biochain.
       # We split the fname by '_' and get the second element.
        ref          = lambda wildcards: config["ref"]["bismark-hisat2"][wildcards.fname.split('_')[1]],
        extra_params = (config["bismark"]["align_hisat2"]["extra_params"]
                        if config["bismark"]["align_hisat2"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        bismark \
            -p {threads} --bam --hisat2 \
            --output_dir result/{wildcards.fname}/{wildcards.trimmer}/bismark-hisat2 \
            --genome_folder {params.ref} \
            -1 {input.r1} -2 {input.r2} {params.extra_params}

        cd result/{wildcards.fname}/{wildcards.trimmer}/bismark-hisat2

        LIB=$(echo "{wildcards.fname}" | cut -d'_' -f1)
        PLATFORM=$(echo "{wildcards.fname}" | cut -d'_' -f4)
        SAMPLE=$(echo "{wildcards.fname}" | cut -d'_' -f2-3)

        samtools sort \
            {wildcards.fname}.R1_bismark_bt2_pe.bam \
            -m 1G |\
        samtools addreplacerg - \
            -r "@RG\\tID:{wildcards.fname}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}\" \
            --output-fmt bam,level=9 \
            -@ {threads} -o {wildcards.fname}.bam

        samtools index -@ {threads} {wildcards.fname}.bam || echo "suppress non-zero exit"

        mv \
            {wildcards.fname}.R1_bismark_bt2_PE_report.txt \
            {wildcards.fname}.align.report

        rm {wildcards.fname}.R1_bismark_bt2_pe.bam
        """
