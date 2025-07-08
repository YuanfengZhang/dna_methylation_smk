configfile: "config/runtime_config.yaml"
from textwrap import dedent

# ! This is unpractical, as same as aryana-bs. Too slow.

rule reckoner:
    input:
        r1            = "result/{BaseName}/{CalibrateParentDir}/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{CalibrateParentDir}/{BaseName}.R2.fq.gz"
    output:
        r1            = "result/{BaseName}/{CalibrateParentDir}/reckoner/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{CalibrateParentDir}/reckoner/{BaseName}.R2.fq.gz",
        r1_log        = "result/{BaseName}/{CalibrateParentDir}/reckoner/{BaseName}.R1.reckoner.log"
    benchmark:
        "result/{BaseName}/{CalibrateParentDir}/reckoner/{BaseName}.calibrate.benchmark"
    params:
        extra_params  = config["reckoner"]["calibrate"]["extra_params"] or ""
    threads: 32
    shell:
        dedent("""\
        export tmp_dir="result/{wildcards.BaseName}/{wildcards.CalibrateParentDir}/reckoner"
        resources/RECKONER/bin/reckoner \\
          -genome 3200000000 \\
          -read {input.r1} -read {input.r2} \\
          -prefix ${{tmp_dir}} \\
          -threads {threads} {params.extra_params}

        mv ${{tmp_dir}}/{wildcards.BaseName}.R1.corrected.fq.gz {output.r1}
        mv ${{tmp_dir}}/{wildcards.BaseName}.R2.corrected.fq.gz {output.r2}
        mv  ${{tmp_dir}}/{wildcards.BaseName}.R1.log {output.r1_log}
        """)
