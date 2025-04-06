configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
export PKG_CONFIG_PATH=/home/zhangyuanfeng/mambaforge/envs/genomic_tools/lib/pkgconfig
make all -j 64
resources/segemehl-0.3.4/segemehl.x -x BL.idx -d BL.fa -t 64
resources/segemehl-0.3.4/segemehl.x \
    -i /hot_warm_data/ref/quartet/DNA/custom_genome/BL/segemehl/BL.idx \
    -d /hot_warm_data/ref/quartet/DNA/custom_genome/BL/segemehl/BL.fa \
    -q input/BS_HF_1_MX.R1.fq.gz \
    -p input/BS_HF_1_MX.R2.fq.gz \
    -t 64 \
    -F 2
"""


def retrieve_method_param(wildcards):
    if wildcards.BaseName.split('_')[0] in ('PS', 'NP'):
        return ""
    else:
        return "-F 2"


rule segemehl_align:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz"
    output:
        bam          = "result/{BaseName}/{AlignParentDir}/segemehl/{BaseName}.bam",
        bai          = "result/{BaseName}/{AlignParentDir}/segemehl/{BaseName}.bam.bai"
    params:
        ref          = lambda wildcards: config["ref"]["segemehl"][wildcards.BaseName.split('_')[1]],
        method_param = lambda wildcards: retrieve_method_param(wildcards),
        extra_params = config["segemehl"]["align"]["extra_params"] or ""
    threads: 8
    conda:
        "../samtools/conda.yaml"
    shell:
        dedent("""
        LIB=$(echo "{wildcards.BaseName}" | cut -d'_' -f1)
        PLATFORM=$(echo "{wildcards.BaseName}" | cut -d'_' -f4)
        SAMPLE=$(echo "{wildcards.BaseName}" | cut -d'_' -f2-3)

        resources/segemehl-0.3.4/segemehl.x \
            {params.ref} -q {input.r1} -p {input.r2} \
            -t {threads} {params.method_param} {params.extra_params} |\
        mbuffer -m 4G -q |\
        samtools sort -@ {threads} - |\
        samtools addreplacerg \
            -r "@RG\\tID:{wildcards.BaseName}\\tSM:${{SAMPLE}}\\tPL:${{PLATFORM}}\\tLB:${{LIB}}" \
            --output-fmt bam,level=9 \
            -@ {threads} -o {output.bam} -

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
