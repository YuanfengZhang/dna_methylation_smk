configfile: "config/runtime_config.yaml"
from textwrap import dedent
"""
abismalidx ./BL.fa BL.abismal_index
export LD_LIBRARY_PATH=/home/zhangyuanfeng/mambaforge/envs/genomic_tools/lib:$LD_LIBRARY_PATH
abismal \
    -i /hot_warm_data/ref/quartet/DNA/custom_genome/BL/abismal/BL.abismal_index \
    -o /mnt/eqa/zhangyuanfeng/methylation/dna_methylation_smk/test/align/abismal.bam \
    -t 16 -B \
    /mnt/eqa/zhangyuanfeng/methylation/dna_methylation_smk/input/BS_HF_1_MX.R1.fq.gz \
    /mnt/eqa/zhangyuanfeng/methylation/dna_methylation_smk/input/BS_HF_1_MX.R2.fq.gz
"""

rule abismal_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz",
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/abismal/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/abismal/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/abismal/{BaseName}.align.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["abismal"][wildcards.BaseName.split('_')[1]],
        extra_params = config["abismal"]["align"]["extra_params"] or ""
    threads: 64
    conda:
        "../samtools/conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d'_' -f1)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d'_' -f4)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        CONDA_LIB=$(which samtools | sed 's/bin\/samtools/lib/')
        export LD_LIBRARY_PATH=$CONDA_LIB:$LD_LIBRARY_PATH
        export LD_PRELOAD=$CONDA_LIB/libhts.so.3

        resources/abismal/bin/abismal \
            -i {params.ref} \
            -o {output.bam}.tmp \
            -t 16 -B {params.extra_params} \
            {input.r1} {input.r2}

        samtools sort -@ {threads} {output.bam}.tmp |\
        samtools addreplacerg \
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \
            --output-fmt bam,level=9 \
            -@ {threads} -o {output.bam} -

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"

        rm {output.bam}.tmp""")
