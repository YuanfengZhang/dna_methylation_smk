configfile: "config/runtime_config.yaml"

rule gatk_bqsr:
    input:
        "result/{BaseName}/{CalibrateParentDir}/{BaseName}.bam"
    output:
        bam          = "result/{BaseName}/{CalibrateParentDir}/gatk-cali/{BaseName}.bam",
        bai          = "result/{BaseName}/{CalibrateParentDir}/gatk-cali/{BaseName}.bam.bai",
        before_table = "result/{BaseName}/{CalibrateParentDir}/gatk-cali/{BaseName}.before.table",
        after_table  = "result/{BaseName}/{CalibrateParentDir}/gatk-cali/{BaseName}.after.table",
        pdf          = "result/{BaseName}/{CalibrateParentDir}/gatk-cali/{BaseName}.AnalyzeCovariates.pdf"
    params:
        snp_vcf      = config["gatk"]["cal_bqsr"]["snp_vcf"],
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        bqsr_params  = (config["gatk"]["cal_bqsr"]["extra_params"]
                        if config["gatk"]["cal_bqsr"]["extra_params"] else ""),
        apply_params = (config["gatk"]["apply_bqsr"]["extra_params"]
                        if config["gatk"]["apply_bqsr"]["extra_params"] else ""),
        cov_params   = (config["gatk"]["analyze_covariates"]["extra_params"]
                        if config["gatk"]["analyze_covariates"]["extra_params"] else ""),
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        gatk \
            --java-options "-Xmx20g -XX:ParallelGCThreads={threads}"\
            BaseRecalibrator \
            -I {input} -O {output.before_table} \
            -R {params.ref} --known-sites {params.snp_vcf} \
            {params.bqsr_params}

        gatk \
            --java-options "-Xmx20g -XX:ParallelGCThreads={threads}"\
            ApplyBQSR \
            -I {input} -O {output.bam} \
            -R {params.ref} \
            --bqsr-recal-file {output.before_table} \
            {params.apply_params}
        
        gatk \
            --java-options "-Xmx20g -XX:ParallelGCThreads={threads}"\
            BaseRecalibrator \
            -I {output.bam} -O {output.after_table} \
            -R {params.ref} --known-sites {params.snp_vcf} \
            {params.bqsr_params}
        
        gatk \
            --java-options "-Xmx20g -XX:ParallelGCThreads={threads}" \
            AnalyzeCovariates \
            -before {output.before_table} \
            -after {output.after_table} \
            -csv {output.pdf} \
            {params.cov_params}
        
        samtools index  -@ {threads} {output.bam} || echo "suppress non-zero exit"
        """
