configfile: "config/runtime_config.yaml"

rule bam2nuc:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam2nuc.txt"
    params:
        ref = lambda wildcards: config["ref"]["bismark-bowtie2"][wildcards.fname.split('_')[1]],
        extra_params = (config["bismark"]["bam2nuc"]["extra_params"]
                        if config["bismark"]["bam2nuc"]["extra_params"] else ""),
        pattern = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 1
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        bam2nuc \
            {wildcards.fname}{params.pattern}.bam \
            --dir . \
            --genome_folder {params.ref} \
            {params.extra_params}
        mv \
            {wildcards.fname}{params.pattern}.nucleotide_stats.txt \
            {wildcards.fname}{params.pattern}.bam2nuc.txt
        """
