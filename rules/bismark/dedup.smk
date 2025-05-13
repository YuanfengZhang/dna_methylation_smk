configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule bismark_dedup:
    input:
        "result/{BaseName}/{DedupParentDir}/{BaseName}.bam"
    output:
        multiext("result/{BaseName}/{DedupParentDir}/bismark-dedup/{BaseName}",
                 ".bam",
                 ".bam.bai")
    benchmark:
        "result/{BaseName}/{DedupParentDir}/bismark-dedup/{BaseName}.dedup.benchmark"
    params:
        extra_params = config["bismark"]["dedup"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.DedupParentDir}
        samtools \\
            sort -n \\
            -@ {threads} \\
            -o bismark-dedup/{wildcards.BaseName}.sort_n_tmp.bam \\
            {wildcards.BaseName}.bam
        cd bismark-dedup
        deduplicate_bismark \\
            --paired --bam \\
            --outfile {wildcards.BaseName} \\
            --output_dir . \\
            {wildcards.BaseName}.sort_n_tmp.bam \\
            {params.extra_params}
        
        samtools sort -@ {threads} \\
            -O bam,level=9 \\
            -o {wildcards.BaseName}.bam \\
            {wildcards.BaseName}.deduplicated.bam

        samtools index -@ {threads} {wildcards.BaseName}.bam || echo "suppress non-zero exit"

        rm {wildcards.BaseName}.sort_n_tmp.bam
        rm {wildcards.BaseName}.deduplicated.bam
        """)
