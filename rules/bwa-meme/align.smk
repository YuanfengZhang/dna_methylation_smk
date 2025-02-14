configfile: "config/runtime_config.yaml"

rule bwa_meme:
    input:
        r1           = "result/{trimmer}/{fname}.R1.fq.gz",
        r2           = "result/{trimmer}/{fname}.R2.fq.gz"
    output:
        bam          = "result/{trimmer}/bwa-meme/{fname}.bam"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-meme"][wildcards.fname.split('_')[1]],
        extra_params = (config["bwa-meme"]["extra_params"]
                        if config["bwa-meme"]["extra_params"] else "")
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        bwa-meme mem \
            -7 -Y -K 100000000 -t {threads} \
            {params.ref} {input.r1} {input.r2} |\
        mbuffer -m 1G -q |\
        samtools sort \
            --output-fmt bam,level=1 \
            -@ {threads} - >{output.bam}
        
        samtools index -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """
