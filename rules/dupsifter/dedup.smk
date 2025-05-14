configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule dupsifter_dedup:
    input:
        bam          = "result/{BaseName}/{DedupParentDir}/{BaseName}.bam"
    output:
        temp("result/{BaseName}/{DedupParentDir}/sort_n_tmp.bam"),
        bam          = "result/{BaseName}/{DedupParentDir}/dupsifter/{BaseName}.bam",
        bai          = "result/{BaseName}/{DedupParentDir}/dupsifter/{BaseName}.bam.bai",
        stats        = "result/{BaseName}/{DedupParentDir}/dupsifter/{BaseName}.dupsifter.stat"
    benchmark:
        "result/{BaseName}/{DedupParentDir}/dupsifter/{BaseName}.dedup.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["dupsifter"][wildcards.BaseName.split('_')[1]],
        extra_params = config["dupsifter"]["dedup"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        tmp_bam="result/{wildcards.BaseName}/{wildcards.DedupParentDir}/sort_n_tmp.bam"
        if [ ! -f $tmp_bam ]; then
            samtools sort -n \\
            -@ {threads} \\
            -o $tmp_bam \\
            {input.bam}
        fi
        dupsifter \\
            -O {output.stats} {params.extra_params} \\
            {params.ref} \\
            $tmp_bam |\\
        mbuffer -m 4G -q |\\
        samtools sort -@ {threads} \\
            -O bam,level=9 \\
            -o {output.bam} -

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
