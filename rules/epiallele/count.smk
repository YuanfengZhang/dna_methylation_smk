configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule epiallele_r:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        copied_bam       = temp("result/{BaseName}/{CountParentDir}/epiallele/{BaseName}.bam"),
        preprocessed_bam = temp("result/{BaseName}/{CountParentDir}/epiallele/{BaseName}.preprocessed.bam"),
        name_sorted_bam  = temp("result/{BaseName}/{CountParentDir}/epiallele/{BaseName}.sortn.bam"),
        raw_report       = temp("result/{BaseName}/{CountParentDir}/epiallele/{BaseName}.raw_report"),
        filtered_report  = "result/{BaseName}/{CountParentDir}/epiallele/{BaseName}.report.zst"
    benchmark:
        "result/{BaseName}/{CountParentDir}/epiallele/{BaseName}.count.benchmark"
    params:
        container    = config["epiallele"]["container"],
        docker_sig   = config["epiallele"]["docker_sig"],
        method       = lambda wildcards: wildcards.BaseName.split('_')[0][: 2],
        ref_dir      = lambda wildcards: config["ref"]["epiallele"]["ref_dir"][wildcards.BaseName.split('_')[1]],
        ref_fname    = lambda wildcards: config["ref"]["epiallele"]["ref_name"][wildcards.BaseName.split('_')[1]],
        extra_params = config["epiallele"]["extra_params"] or ""
    threads: 32
    conda:
        "conda.yaml"
    shell:
        dedent("""
        export tmp_dir="$(pwd)/result/{wildcards.BaseName}/{wildcards.CountParentDir}"
        mkdir -p ${{tmp_dir}}/epiallele

        tmp_read="$(samtools view {input} | head -n 1 | awk '{{gsub("\\t", "_"); print}}')" || echo ""

        python rules/epiallele/prepare.py \\
          -i ${{tmp_dir}} \\
          -b {wildcards.BaseName} \\
          -r "$tmp_read" \\
          -t {threads}

        if [ {params.docker_sig} == "docker" ]; then
          mount_flag="-v"
          container_run="docker run --rm"
        else
          mount_flag="-B"
          container_run="singularity exec"
        fi

        ${{container_run}} \\
          ${{mount_flag}} {params.ref_dir}:/data/ref \\
          ${{mount_flag}} ${{tmp_dir}}/epiallele:/data/process \\
          ${{mount_flag}} $(pwd)/rules/epiallele:/data/script \\
          --cpus {threads} \\
          {params.container} \\
          Rscript /data/script/epiallele.R \\
            --basename {wildcards.BaseName} \\
            --ref-name {params.ref_fname} \\
            --threads {threads} \\
            {params.extra_params}

        awk -F'\t' '
           NR == 1 || (NR > 1 && ($5 + $6) >= 5)
          ' {output.raw_report} |\\
        zstd -T{threads} -19 > {output.filtered_report}
        """)
