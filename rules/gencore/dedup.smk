configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule gencore_dedup:
    input:
        bam              = "result/{BaseName}/{DedupParentDir}/{BaseName}.bam"
    output:
        bam              = "result/{BaseName}/{DedupParentDir}/gencore/{BaseName}.bam",
        bai              = "result/{BaseName}/{DedupParentDir}/gencore/{BaseName}.bam.bai",
        json             = "result/{BaseName}/{DedupParentDir}/gencore/{BaseName}.gencore.json",
        html             = "result/{BaseName}/{DedupParentDir}/gencore/{BaseName}.gencore.html"
    benchmark:
        "result/{BaseName}/{DedupParentDir}/gencore/{BaseName}.dedup.benchmark"
    params:
        ref              = config["ref"]["gencore"][wildcards.BaseName.split('_')[1]],
        supporting_reads = config["gencore"]["dedup"]["supporting_reads"],
        score_threshold  = config["gencore"]["dedup"]["score_threshold"], 
        ratio_threshold  = config["gencore"]["dedup"]["ratio_threshold"], 
        extra_params     = lambda wildcards: config["gencore"]["dedup"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        gencore -i {input.bam} \\
            -r {params.ref} \\
            --supporting_reads {params.supporting_reads} \\
            --score_threshold {params.score_threshold} \\
            --ratio_threshold {params.ratio_threshold} \\
            -j {output.json} \\
            -h {output.html} {params.extra_params} |\\
        mbuffer -m 4G -q |\\
        samtools sort -@ {threads} \\
            -O bam,level=9 \\
            -o {output.bam} -

        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """)

