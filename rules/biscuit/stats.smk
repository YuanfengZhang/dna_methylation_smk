configfile: "config/runtime_config.yaml"

rule biscuit_qc:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        multiext("result/{fname}/{trimmer}/{aligner}/{deduper}/"
                 "{fname}{bqsr}",
                 "biscuit.CpGRetentionByReadPos.txt",
                 "biscuit.CpHRetentionByReadPos.txt",
                 "biscuit.dup_report.txt",
                 "biscuit.isize_table.txt",
                 "biscuit.mapq_table.txt",
                 "biscuit.strand_table.txt",
                 "biscuit.totalReadConversionRate.txt")
    params:
        ref          = lambda wildcards: config["ref"]["biscuit"][wildcards.fname.split('_')[1]],
        qc_params    = (config["biscuit"]["qc"]["extra_params"]
                        if config["biscuit"]["qc"]["extra_params"] else ""),
        pattern      = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 1
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        biscuit qc \
            {params.ref} \
            {wildcards.fname}{params.pattern}.bam \
            {wildcards.fname}{params.pattern} {params.qc_params}
        
        mv \
            {wildcards.fname}{params.pattern}_CpGRetentionByReadPos.txt\
            {wildcards.fname}{params.pattern}.biscuit.CpGRetentionByReadPos.txt
        mv \
            {wildcards.fname}{params.pattern}_CpHRetentionByReadPos.txt\
            {wildcards.fname}{params.pattern}.biscuit.CpHRetentionByReadPos.txt
        mv \
            {wildcards.fname}{params.pattern}_dup_report.txt\
            {wildcards.fname}{params.pattern}.biscuit.dup_report.txt
        mv \
            {wildcards.fname}{params.pattern}_isize_table.txt\
            {wildcards.fname}{params.pattern}.biscuit.isize_table.txt
        mv \
            {wildcards.fname}{params.pattern}_mapq_table.txt\
            {wildcards.fname}{params.pattern}.biscuit.mapq_table.txt
        mv \
            {wildcards.fname}{params.pattern}_strand_table.txt\
            {wildcards.fname}{params.pattern}.biscuit.strand_table.txt
        mv \
            {wildcards.fname}{params.pattern}_totalReadConversionRate.txt\
            {wildcards.fname}{params.pattern}.biscuit.totalReadConversionRate.txt
        """

rule biscuit_bsstrand:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        ("result/{fname}/{trimmer}/{aligner}/{deduper}/"
         "{fname}{bqsr}biscuit.bsstrand.txt")
    params:
        ref             = lambda wildcards: config["ref"]["biscuit"][wildcards.fname.split('_')[1]],
        bsstrand_params = (config["biscuit"]["bsstrand"]["extra_params"]
                           if config["biscuit"]["bsstrand"]["extra_params"] else ""),
        pattern         = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 1
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        biscuit bsstrand \
            {params.ref} {wildcards.fname}{params.pattern}.bam {params.bsstrand_params} > {wildcards.fname}{params.pattern}.biscuit.bsstrand.txt 2>&1
        """

rule biscuit_cinread:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        ("result/{fname}/{trimmer}/{aligner}/{deduper}/"
         "{fname}{bqsr}biscuit.cinread.txt")
    params:
        ref            = lambda wildcards: config["ref"]["biscuit"][wildcards.fname.split('_')[1]],
        cinread_params = (config["biscuit"]["cinread"]["extra_params"]
                          if config["biscuit"]["cinread"]["extra_params"] else ""),
        pattern        = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 1
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        biscuit cinread \
            {params.ref} {wildcards.fname}{params.pattern}.bam {params.cinread_params}\
            -o {wildcards.fname}{params.pattern}.biscuit.cinread.txt
        """
