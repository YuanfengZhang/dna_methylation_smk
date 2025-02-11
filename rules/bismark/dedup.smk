configfile: "config/runtime_config.yaml"

rule bismark_dedup:
    input:
        "result/{fname}/{trimmer}/{aligner}/{fname}.bam"
    output:
        multiext("result/{fname}/{trimmer}/{aligner}/bismark/{fname}",
                 ".bam",
                 ".bam.bai")
    params:
        extra_params = (config["bismark"]["dedup"]["extra_params"]
                        if config["bismark"]["dedup"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}
        samtools \
            sort -n \
            -@ {threads} \
            -o bismark/{wildcards.fname}.sort_n_tmp.bam \
            {wildcards.fname}.bam
        cd bismark
        deduplicate_bismark \
            --paired --bam \
            --outfile {wildcards.fname} \
            --output_dir . \
            {wildcards.fname}.sort_n_tmp.bam \
            {params.extra_params}
        
        samtools sort -@ {threads} \
            -O bam,level=9 \
            -o {wildcards.fname}.bam \
            {wildcards.fname}.deduplicated.bam

        samtools index -@ {threads} {wildcards.fname}.bam

        rm {wildcards.fname}.sort_n_tmp.bam
        rm {wildcards.fname}.deduplicated.bam
        """
