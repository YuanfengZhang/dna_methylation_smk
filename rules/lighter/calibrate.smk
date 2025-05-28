configfile: "config/runtime_config.yaml"
from textwrap import dedent


def sampling_rate(wildcards):
    lab: str = wildcards['BaseName'].split('_')[0]
    avg_depth_dict: dict[str, int] = {
        'BS1': 30, 'BS2': 30, 'BS3': 30, 'BS4': 30,
        'EM1': 30, 'EM2': 30, 'EM3': 30, 'EM4': 30,
        'PS1': 30, 'PS2': 30, 'PS3': 30, 'RR1': 30}
    
    return round(avg_depth_dict.get(lab, 30), 4)


rule lighter:
    input:
        r1            = "result/{BaseName}/{CalibrateParentDir}/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{CalibrateParentDir}/{BaseName}.R2.fq.gz"
    output:
        r1            = "result/{BaseName}/{CalibrateParentDir}/lighter/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{CalibrateParentDir}/lighter/{BaseName}.R2.fq.gz"
    benchmark:
        "result/{BaseName}/{CalibrateParentDir}/lighter/{BaseName}.calibrate.benchmark"
    params:
        sampling_rate = lambda wildcards: sampling_rate(wildcards),
        extra_params  = config["lighter"]["calibrate"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""\
        export tmp_dir="result/{wildcards.BaseName}/{wildcards.CalibrateParentDir}/lighter"
        lighter \\
          -r {input.r1} -r {input.r2} \\
          -od ${{tmp_dir}} \\
          -t {threads} -k 19 3200000000 {params.sampling_rate} {params.extra_params}

        mv ${{tmp_dir}}/{wildcards.BaseName}.R1.cor.fq.gz {output.r1}
        mv ${{tmp_dir}}/{wildcards.BaseName}.R2.cor.fq.gz {output.r2}
        """)
