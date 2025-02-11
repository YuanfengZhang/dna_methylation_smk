configfile: "config/runtime_config.yaml"

rule fake_dedup:
    input:
        bam = "result/{fname}/{trimmer}/{aligner}/{fname}.bam",
        bai = "result/{fname}/{trimmer}/{aligner}/{fname}.bam.bai"
    output:
        bam = "result/{fname}/{trimmer}/{aligner}/no_dedup/{fname}.bam",
        bai = "result/{fname}/{trimmer}/{aligner}/no_dedup/{fname}.bam.bai"
    threads: 1
    shell:
        """
        ln -s "$(pwd)/{input.bam}" "$(pwd)/{output.bam}"
        ln -s "$(pwd)/{input.bai}" "$(pwd)/{output.bai}"
        """
