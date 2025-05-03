configfile: "config/runtime_config.yaml"
from textwrap import dedent
"""
pycoMeth CGI_Finder -f BL.fa -b BL.cgi.bed
prepare_genomes.sh BL.fa BL.cgi.bed ./aryana .

!!! The builtin prepare_genomes.sh is troublesome.
!!! It exits before finishing the bwa index, leading to truncated index files on cluster.

aryana_bs \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/aryana/BL.fa \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/aryana/aryana \
    /hot_warm_data/ref/quartet/DNA/custom_genome/BL/aryana/BL.cgi.bed \
    merged.fq.gz \
    /mnt/eqa/zhangyuanfeng/methylation/dna_methylation_smk/test/align/aryana \
    ar="-p 12"
"""

rule aryana_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/aryana/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/aryana/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/aryana/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["aryana"][wildcards.BaseName.split('_')[1]],
        extra_params = config["aryana"]["align"]["extra_params"] or ""
    threads: 64
    conda:
        "../samtools/conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d _ -f1 | cut -c 1-2)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d _ -f1)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        parent_dir="result/{wildcards.BaseName}/{wildcards.AlignParentDir}"

        seqtk mergepe \\
            ${{parent_dir}}/{wildcards.BaseName}.R1.fq.gz \\
            ${{parent_dir}}/{wildcards.BaseName}.R2.fq.gz \\
            >${{parent_dir}}/{wildcards.BaseName}.fq

        resources/aryana/aryana_bs \\
            {params.ref} \\
            ${{parent_dir}}/{wildcards.BaseName}.fq \\
            ${{parent_dir}}/aryana/{wildcards.BaseName} \\
            ar='-p {threads}'

        samtools sort -@ {threads} ${{parent_dir}}/aryana/{wildcards.BaseName}.sam |\\
        samtools addreplacerg \\
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \\
            --output-fmt bam,level=9 \\
            -@ {threads} -o  ${{parent_dir}}/aryana/{wildcards.BaseName}.bam -

        samtools index \\
            -@ {threads} \\
            ${{parent_dir}}/aryana/{wildcards.BaseName}.bam || echo "suppress non-zero exit"

        rm ${{parent_dir}}/{wildcards.BaseName}.fq ${{parent_dir}}/aryana/{wildcards.BaseName}.sam
        """)
