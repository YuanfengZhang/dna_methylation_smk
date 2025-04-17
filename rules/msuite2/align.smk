configfile: "config/runtime_config.yaml"
from textwrap import dedent


def get_ref_label(BaseName) -> str:
    label_map = {'D5': 'D5', 'D6': 'D5',
                 'F7': 'F7', 'M8': 'M8',
                 'BC': 'BC', 'BL': 'BL',
                 'T1': 'D5', 'T2': 'D5',
                 'T3': 'D5', 'T4': 'D5',
                 'HF': 'BL'}
    return label_map.get(BaseName.split('_')[1], 'BL')


def get_method_param(BaseName) -> str:
    method_map = {'BS': '-3 -m BS',
                  'EM': '-3 -m BS',
                  'RR': '-3 -m BS',
                  'TB': '-3 -m BS',
                  'PS': '-4 -m TAPS'}
    return method_map.get(BaseName.split('_')[0][: 2], '-3 -m BS')


def get_platform_param(BaseName) -> str:
    platform_map = {'BS1': 'bgi', 'BS2': 'illumina', 'BS3': 'illumina', 'BS4': 'illumina',
                    'EM1': 'illumina', 'EM2': 'bgi', 'EM3': 'bgi', 'EM4': 'bgi',
                    'PS1': 'bgi', 'PS2': 'bgi', 'PS3': 'bgi',
                    'RR1': 'bgi', 'TB1': 'bgi'}
    return platform_map.get(BaseName.split('_')[0], 'illumina')


rule msuite2_bowtie2:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz"
    output:
        "result/{BaseName}/{AlignParentDir}/msuite2-bowtie2/{BaseName}.bam",
        "result/{BaseName}/{AlignParentDir}/msuite2-bowtie2/{BaseName}.bam.bai",
        "result/{BaseName}/{AlignParentDir}/msuite2-bowtie2/{BaseName}.msuite2.bedgraph.gz",
        "result/{BaseName}/{AlignParentDir}/msuite2-bowtie2/Msuite2.report/index.html"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/msuite2-bowtie2/{BaseName}.align.benchmark"
    params:
        label        = lambda wildcards: get_ref_label(wildcards.BaseName),
        platform     = lambda wildcards: get_platform_param(wildcards.BaseName),
        method_param = lambda wildcards: get_method_param(wildcards.BaseName),
        ref          = lambda wildcards: config["ref"]["msuite2"][wildcards.BaseName.split('_')[1]],
        extra_params = config["msuite2"]["extra_params"]
                       if config["msuite2"]["extra_params"] else ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        output_dir="result/{wildcards.BaseName}/{wildcards.AlignParentDir}/msuite2-bowtie2"

        [ -L "resources/Msuite2/index/{params.label}" ] || ln -s {params.ref} "$(pwd)/resources/Msuite2/index/{params.label}"
        rm -rf result/{wildcards.BaseName}/{wildcards.AlignParentDir}/msuite2-bowtie2
        resources/Msuite2/msuite2 \
            -x {params.label} -k {params.platform} \
            {params.method_param} \
            -1 {input.r1} -2 {input.r2} \
            -o ${{output_dir}} \
            --aligner bowtie2 \
            -p {threads} {params.extra_params}
        make -C "${{output_dir}}"

        bash rules/msuite2/format.sh ${{output_dir}} {wildcards.BaseName} {threads}
        """)

rule msuite2_hisat2:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz"
    output:
        "result/{BaseName}/{AlignParentDir}/msuite2-hisat2/{BaseName}.bam",
        "result/{BaseName}/{AlignParentDir}/msuite2-hisat2/{BaseName}.bam.bai",
        "result/{BaseName}/{AlignParentDir}/msuite2-hisat2/{BaseName}.msuite2.bedgraph.gz",
        "result/{BaseName}/{AlignParentDir}/msuite2-hisat2/Msuite2.report/index.html"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/msuite2-bowtie2/{BaseName}.align.benchmark"
    params:
        label        = lambda wildcards: get_ref_label(wildcards.BaseName),
        platform     = lambda wildcards: get_platform_param(wildcards.BaseName),
        method_param = lambda wildcards: get_method_param(wildcards.BaseName),
        ref          = lambda wildcards: config["ref"]["msuite2"][wildcards.BaseName.split('_')[1]],
        extra_params = config["msuite2"]["extra_params"]
                       if config["msuite2"]["extra_params"] else ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        output_dir="result/{wildcards.BaseName}/{wildcards.AlignParentDir}/msuite2-hisat2"

        [ -L "resources/Msuite2/index/{params.label}" ] || ln -s {params.ref} "$(pwd)/resources/Msuite2/index/{params.label}"
        rm -rf result/{wildcards.BaseName}/{wildcards.AlignParentDir}/msuite2-hisat2
        resources/Msuite2/msuite2 \
            -x {params.label} -k {params.platform} \
            {params.method_param} \
            -1 {input.r1} -2 {input.r2} \
            -o ${{output_dir}} \
            --aligner hisat2 \
            -p {threads} {params.extra_params}
        make -C ${{output_dir}}

        bash rules/msuite2/format.sh ${{output_dir}} {wildcards.BaseName} {threads}
        """)
