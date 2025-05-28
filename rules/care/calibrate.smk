configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule care:
    input:
        r1            = "result/{BaseName}/{CalibrateParentDir}/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{CalibrateParentDir}/{BaseName}.R2.fq.gz"
    output:
        tmp_r1        = temp("result/{BaseName}/{CalibrateParentDir}/care/{BaseName}.calibrated.R1.fq"),
        tmp_r2        = temp("result/{BaseName}/{CalibrateParentDir}/care/{BaseName}.calibrated.R2.fq"),
        r1            = "result/{BaseName}/{CalibrateParentDir}/care/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{CalibrateParentDir}/care/{BaseName}.R2.fq.gz"
    benchmark:
        "result/{BaseName}/{CalibrateParentDir}/care/{BaseName}.calibrate.benchmark"
    params:
        device        = config["care"]["calibrate"]["device"],  # cpu / gpu
        max_mem       = config["care"]["calibrate"]["max_mem"],  # 48
        avg_cov       = config["care"]["calibrate"]["avg_cov"],  # 30
        extra_params  = config["care"]["calibrate"]["extra_params"] or ""
    threads: 8
    conda:
        "../biscuit/conda.yaml"
    shell:
        dedent("""\
        tmp_dir="result/{wildcards.BaseName}/{wildcards.CalibrateParentDir}/care"
        mkdir -p ${{tmp_dir}}
        ../resources/CARE/care-{params.device} \\
          -i {input.r1} --pairmode PE \\
          -d ${{tmp_dir}} -o {wildcards.BaseName}.calibrated.R1.fq \\
          -c {params.avg_cov} \\
          -m {params.max_mem}G -t {threads} \\
          -p {params.extra_params}

        ../resources/CARE/care-{params.device} \\
          -i {input.r2} --pairmode PE \\
          -d ${{tmp_dir}} -o {wildcards.BaseName}.calibrated.R2.fq \\
          -c {params.avg_cov} \\
          -m {params.max_mem}G -t {threads} \\
          -p {params.extra_params}
        
        pigz --best -c -p {threads} {output.tmp_r1} > {output.r1}
        pigz --best -c -p {threads} {output.tmp_r2} > {output.r2}
        """)
