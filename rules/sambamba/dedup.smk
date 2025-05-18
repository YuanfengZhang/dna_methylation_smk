configfile: "config/runtime_config.yaml"
from textwrap import dedent


# memory consumption grows by 2Gb per each 100M reads

rule sambamba_dedup:
    input:
        bam                = "result/{BaseName}/{DedupParentDir}/{BaseName}.bam"
    output:
        tmp_bam            = temp("result/{BaseName}/{DedupParentDir}/sambamba/tmp.bam"),
        bam                = "result/{BaseName}/{DedupParentDir}/sambamba/{BaseName}.bam",
        bai                = "result/{BaseName}/{DedupParentDir}/sambamba/{BaseName}.bam.bai"
    benchmark:
        "result/{BaseName}/{DedupParentDir}/sambamba/{BaseName}.dedup.benchmark"
    params:
        ref                = lambda wildcards: config["ref"]["sambamba"][wildcards.BaseName.split('_')[1]],
        hash_table_size    = lambda wildcards: config["sambamba"]["dedup"]["hash_table_size"],
        overflow_list_size = lambda wildcards: config["sambamba"]["dedup"]["overflow_list_size"],
        sort_buffer_size   = lambda wildcards: config["sambamba"]["dedup"]["sort_buffer_size"],
        extra_params       = lambda wildcards: config["sambamba"]["dedup"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        sambamba markdup \\
            --hash-table-size {params.hash_table_size} \\
            --overflow-list-size {params.overflow_list_size} \\
            --sort-buffer-size {params.sort_buffer_size} \\
            -t {threads} {params.extra_params} \\
            -l 0 {input.bam} {output.tmp_bam} 
        samtools sort -@ {threads} \\
            -O bam,level=9 \\
            -o {output.bam} \\
            {output.tmp_bam}

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)
