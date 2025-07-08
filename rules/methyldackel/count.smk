configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
methyldackel output:
track type="bedGraph" description="methyldackel/PS3_D6_1 CpG methylation levels"
chr1    26932   26933   0       0       1
chr1    26957   26958   0       0       1

with header.
chrom, start, end, beta, m, u.

For TAPS, we need to flip the m and u columns.
"""

rule methyldackel_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        tmp_c        = temp("result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}_CpG.bedGraph"),
        c_bedgraph   = "result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.bedgraph.gz",
        tmp_cpg      = temp("result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.bedGraph"),
        cpg_bedgraph = "result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.merged.bedgraph.gz",
    benchmark:
        "result/{BaseName}/{CountParentDir}/methyldackel/{BaseName}.count.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        method       = lambda wildcards: wildcards.BaseName.split('_')[0][: 2],
        extra_params = config["methyldackel"]["count"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        export tmp_dir="result/{wildcards.BaseName}/{wildcards.CountParentDir}"
        mkdir -p ${{tmp_dir}}/methyldackel
        MethylDackel extract \\
            {params.ref} {input} \\
            -o ${{tmp_dir}}/methyldackel/{wildcards.BaseName} \\
            -@ {threads} {params.extra_params}

        MethylDackel mergeContext \\
            {params.ref} {output.tmp_c} \\
            -o {output.tmp_cpg} {params.extra_params}

        if [ {params.method} = "PS" ]; then
          awk 'BEGIN{{FS=OFS="\\t"}} NR > 1 && ($5 + $6) >= 5 {{
            print $1, $2, $3, $6, $5}}' {output.tmp_c} |\
          pigz -p {threads} --best > {output.c_bedgraph}
          awk 'BEGIN{{FS=OFS="\\t"}} NR > 1 && ($5 + $6) >= 5 {{
            print $1, $2, $3, $6, $5}}' {output.tmp_cpg} |\
          pigz -p {threads} --best > {output.cpg_bedgraph}
        else
          awk 'BEGIN{{FS=OFS="\\t"}} NR > 1 && ($5 + $6) >= 5 {{
            print $1, $2, $3, $5, $6}}' {output.tmp_c} |\
          pigz -p {threads} --best > {output.c_bedgraph}

          awk 'BEGIN{{FS=OFS="\\t"}} NR > 1 && ($5 + $6) >= 5 {{
            print $1, $2, $3, $5, $6}}' {output.tmp_cpg} |\
          pigz -p {threads} --best > {output.cpg_bedgraph}
        fi
        """)
