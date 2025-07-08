configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule care:
    input:
        r1            = "result/{BaseName}/{CalibrateParentDir}/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{CalibrateParentDir}/{BaseName}.R2.fq.gz"
    output:
        tmp_output    = temp("result/{BaseName}/{CalibrateParentDir}/care/{BaseName}.calibrated.fq"),
        r1            = "result/{BaseName}/{CalibrateParentDir}/care/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{CalibrateParentDir}/care/{BaseName}.R2.fq.gz"
    benchmark:
        "result/{BaseName}/{CalibrateParentDir}/care/{BaseName}.calibrate.benchmark"
    params:
        device        = config["care"]["calibrate"]["device"],  # cpu / gpu
        max_mem       = config["care"]["calibrate"]["max_mem"],  # 48
        avg_cov       = config["care"]["calibrate"]["avg_cov"],  # 30
        extra_params  = config["care"]["calibrate"]["extra_params"] or ""
    threads: 32
    conda:
        "../trie-dedup/conda.yaml"
    shell:
        dedent("""\
        tmp_dir="result/{wildcards.BaseName}/{wildcards.CalibrateParentDir}/care"
        mkdir -p ${{tmp_dir}}

        resources/CARE/care-{params.device} \\
          -i {input.r1} -i {input.r2} --pairmode PE \\
          -d ${{tmp_dir}} -o {wildcards.BaseName}.calibrated.fq \\
          -c {params.avg_cov} \\
          -m {params.max_mem}G -t {threads} \\
          -p {params.extra_params}
        
        first_line="$(head -n1 {output.tmp_output})"
        mode1_regex='^@.*/[12]$'
        mode2_regex='^@.* [12]:'

        if [[ "$first_line" =~ $mode1_regex ]]; then
            r1_pattern='^@.*/1$'
            r2_pattern='^@.*/2$'
        elif [[ "$first_line" =~ $mode2_regex ]]; then
            r1_pattern='^@.* 1:'
            r2_pattern='^@.* 2:'
        else
            echo "Error: First line does not match required patterns"
            echo "Must end with /1 or /2, or contain ' 1:' or ' 2:'"
            exit 1
        fi

        grep -A 3 --no-group-separator "$r1_pattern" {output.tmp_output} | \
        pigz --best -p {threads} > {output.r1} || {{ echo "Error processing R1 reads"; exit 1; }}

        grep -A 3 --no-group-separator "$r2_pattern" {output.tmp_output} | \
        pigz --best -p {threads} > {output.r2} || {{ echo "Error processing R2 reads"; exit 1; }}
        """)
