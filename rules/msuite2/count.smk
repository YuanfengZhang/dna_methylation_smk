configfile: "config/runtime_config.yaml"


def get_ref_label(fname) -> str:
    label_map = {'D5': 'D5', 'D6': 'D5',
                 'F7': 'F7', 'M8': 'M8',
                 'BC': 'BC', 'BL': 'BL',
                 'T1': 'D5', 'T2': 'D5',
                 'T3': 'D5', 'T4': 'D5',
                 'HF': 'BL'}
    return label_map.get(fname.split('_')[1], 'BL')


def get_method_param(fname) -> str:
    method_map = {'BS': '-3 -m BS',
                  'EM': '-3 -m BS',
                  'RR': '-3 -m BS',
                  'TB': '-3 -m BS',
                  'PS': '-4 -m TAPS'}
    return method_map.get(fname.split('_')[0], '-3 -m BS')


def get_platform_param(fname) -> str:
    platform_map = {'Biochain': 'illumina', 'GENESEEQ': 'illumina',
                    'MGI': 'bgi', 'HaploX': 'illumina',
                    'FDU': 'bgi', 'HIT': 'bgi',
                    'singlera': 'illumina', 'IPM-GBA': 'bgi',
                    'FAHJLU': 'bgi', 'FAHZZU': 'bgi', 'Geneplus': 'bgi',
                    'SnT': 'bgi', 'iGeneTech': 'bgi'}
    return platform_map.get(fname.split('_')[3], 'illumina')


rule msuite2_bowtie2:
    input:
        r1           = "result/{fname}/{trimmer}/{fname}.R1.fq.gz",
        r2           = "result/{fname}/{trimmer}/{fname}.R2.fq.gz"
    output:
        "result/{fname}/{trimmer}/msuite2-bowtie2/{fname}.bam",
        "result/{fname}/{trimmer}/msuite2-bowtie2/{fname}.msuite2.bedgraph.gz",
        "result/{fname}/{trimmer}/msuite2-bowtie2/Msuite2.report/index.html"
    params:
        label        = lambda wildcards: get_ref_label(wildcards.fname),
        platform     = lambda wildcards: get_platform_param(wildcards.fname),
        method_param = lambda wildcards: get_method_param(wildcards.fname),
        ref          = lambda wildcards: config["ref"]["msuite2"][wildcards.fname.split('_')[1]],
        extra_params = (config["msuite2"]["extra_params"]
                        if config["msuite2"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        output_dir="result/{wildcards.fname}/{wildcards.trimmer}/msuite2-bowtie2"

        [ -L "resources/Msuite2/index/{params.label}" ] || ln -s {params.ref} "$(pwd)/resources/Msuite2/index/{params.label}"
        rm -rf result/{wildcards.fname}/{wildcards.trimmer}/msuite2-bowtie2
        resources/Msuite2/msuite2 \
            -x {params.label} -k {params.platform} \
            {params.method_param} \
            -1 {input.r1} -2 {input.r2} \
            -o ${{output_dir}} \
            --aligner bowtie2 \
            -p {threads} {params.extra_params}
        make -C "${{output_dir}}"

        bash rules/msuite2/format.sh ${{output_dir}} {wildcards.fname} {threads}
        """

rule msuite2_hisat2:
    input:
        r1           = "result/{fname}/{trimmer}/{fname}.R1.fq.gz",
        r2           = "result/{fname}/{trimmer}/{fname}.R2.fq.gz"
    output:
        "result/{fname}/{trimmer}/msuite2-hisat2/{fname}.bam",
        "result/{fname}/{trimmer}/msuite2-hisat2/{fname}.msuite2.bedgraph.gz",
        "result/{fname}/{trimmer}/msuite2-hisat2/Msuite2.report/index.html"
    params:
        label        = lambda wildcards: get_ref_label(wildcards.fname),
        platform     = lambda wildcards: get_platform_param(wildcards.fname),
        method_param = lambda wildcards: get_method_param(wildcards.fname),
        ref          = lambda wildcards: config["ref"]["msuite2"][wildcards.fname.split('_')[1]],
        extra_params = (config["msuite2"]["extra_params"]
                        if config["msuite2"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        output_dir="result/{wildcards.fname}/{wildcards.trimmer}/msuite2-hisat2"

        [ -L "resources/Msuite2/index/{params.label}" ] || ln -s {params.ref} "$(pwd)/resources/Msuite2/index/{params.label}"
        rm -rf result/{wildcards.fname}/{wildcards.trimmer}/msuite2-hisat2
        resources/Msuite2/msuite2 \
            -x {params.label} -k {params.platform} \
            {params.method_param} \
            -1 {input.r1} -2 {input.r2} \
            -o ${{output_dir}} \
            --aligner hisat2 \
            -p {threads} {params.extra_params}
        make -C ${{output_dir}}

        bash rules/msuite2/format.sh ${{output_dir}} {wildcards.fname} {threads}
        """
