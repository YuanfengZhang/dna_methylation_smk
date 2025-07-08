configfile: "config/runtime_config.yaml"
from textwrap import dedent


# ! the latest released dnmtools cannot skip the reads with $3 (chromosome) that not in the ref.fa,
# ! See https://github.com/smithlabcode/dnmtools/issues/260
# ! The current solution is to use awk to get the legal records, run dnmtools count, and then remove it.

rule dnmtools_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        temp_sam      = temp("result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.mapped.sam"),
        temp_bedgraph = temp("result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.dnmtools.tmp"),
        bedgraph      = "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.dnmtools.gz"
    benchmark:
        "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.count.benchmark"
    params:
        ref           = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        method        = lambda wildcards: wildcards.BaseName.split('_')[0][: 2],
        extra_params  = config["dnmtools"]["count"]["extra_params"] or ""
    threads: 32
    conda:
        "conda.yaml"
    shell:
        dedent("""
        tmp_dir="result/{wildcards.BaseName}/{wildcards.CountParentDir}"
        mkdir -p "$tmp_dir/dnmtools"

        echo "0" > "$tmp_dir/dnmtools.log"
        if [ ! -f {output.temp_sam} ]; then
          samtools view -h {input} |\\
          awk '$3 != "*"' > {output.temp_sam}.tmp
          dnmtools uniq \\
            -t {threads} \\
            {output.temp_sam}.tmp \\
            {output.temp_sam}
        fi
        echo "1" >> "$tmp_dir/dnmtools/dnmtools.log"
        dnmtools counts \\
          -c {params.ref} -t {threads} \\
          -progress {params.extra_params} \\
          -o {output.temp_bedgraph} \\
          {output.temp_sam}
        echo "2" >> "$tmp_dir/dnmtools/dnmtools.log"
        if [ {params.method} = "PS" ]; then
          mv {output.temp_bedgraph} {output.temp_bedgraph}.tmp
          awk 'BEGIN {{FS=OFS="\\t"}} {{$5 = 1 - $5}} 1' \\
            {output.temp_bedgraph}.tmp > {output.temp_bedgraph}
          rm {output.temp_bedgraph}.tmp
        fi
        echo "3" >> "$tmp_dir/dnmtools/dnmtools.log"
        awk -F'\\t' '$6 > 4' \\
          {output.temp_bedgraph} \\
          | pigz -p {threads} --best > {output.bedgraph}
        """)


rule dnmtools_levels:
    input:
        "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.dnmtools.tmp"
    output:
        "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.dnmtools.stats"
    benchmark:
        "result/{BaseName}/{CountParentDir}/dnmtools/{BaseName}.level.benchmark"
    params:
        extra_params = config["dnmtools"]["levels"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        dnmtools levels -o {output} {input} {params.extra_params}
        """)
