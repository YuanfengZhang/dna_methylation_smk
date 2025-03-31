configfile: "config/runtime_config.yaml"

rule biscuit_qc:
    input:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.bam"
    output:
        multiext("result/{BaseName}/{BamStatsParentDir}/"
                 "{BaseName}",
                 ".biscuit.CpGRetentionByReadPos.txt",
                 ".biscuit.CpHRetentionByReadPos.txt",
                 ".biscuit.dup_report.txt",
                 ".biscuit.isize_table.txt",
                 ".biscuit.mapq_table.txt",
                 ".biscuit.strand_table.txt",
                 ".biscuit.totalReadConversionRate.txt")
    params:
        ref          = lambda wildcards: config["ref"]["biscuit"][wildcards.BaseName.split('_')[1]],
        qc_params    = (config["biscuit"]["qc"]["extra_params"]
                        if config["biscuit"]["qc"]["extra_params"] else "")
    threads: 1
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.BaseName}/{wildcards.BamStatsParentDir}
        biscuit qc \
            {params.ref} \
            {wildcards.BaseName}.bam \
            {wildcards.BaseName} {params.qc_params}
        
        mv \
            {wildcards.BaseName}_CpGRetentionByReadPos.txt\
            {wildcards.BaseName}.biscuit.CpGRetentionByReadPos.txt
        mv \
            {wildcards.BaseName}_CpHRetentionByReadPos.txt\
            {wildcards.BaseName}.biscuit.CpHRetentionByReadPos.txt
        mv \
            {wildcards.BaseName}_dup_report.txt\
            {wildcards.BaseName}.biscuit.dup_report.txt
        mv \
            {wildcards.BaseName}_isize_table.txt\
            {wildcards.BaseName}.biscuit.isize_table.txt
        mv \
            {wildcards.BaseName}_mapq_table.txt\
            {wildcards.BaseName}.biscuit.mapq_table.txt
        mv \
            {wildcards.BaseName}_strand_table.txt\
            {wildcards.BaseName}.biscuit.strand_table.txt
        mv \
            {wildcards.BaseName}_totalReadConversionRate.txt\
            {wildcards.BaseName}.biscuit.totalReadConversionRate.txt
        """

rule biscuit_bsstrand:
    input:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.bam"
    output:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.biscuit.bsstrand.txt"
    params:
        ref             = lambda wildcards: config["ref"]["biscuit"][wildcards.BaseName.split('_')[1]],
        bsstrand_params = (config["biscuit"]["bsstrand"]["extra_params"]
                           if config["biscuit"]["bsstrand"]["extra_params"] else "")
    threads: 1
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.BaseName}/{wildcards.BamStatsParentDir}
        biscuit bsstrand \
            {params.ref} {wildcards.BaseName}.bam {params.bsstrand_params} > {wildcards.BaseName}.biscuit.bsstrand.txt 2>&1
        """

rule biscuit_cinread:
    input:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.bam"
    output:
        ("result/{BaseName}/{BamStatsParentDir}/"
         "{BaseName}.biscuit.cinread.txt")
    params:
        ref            = lambda wildcards: config["ref"]["biscuit"][wildcards.BaseName.split('_')[1]],
        cinread_params = (config["biscuit"]["cinread"]["extra_params"]
                          if config["biscuit"]["cinread"]["extra_params"] else "")
    threads: 1
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.BaseName}/{wildcards.BamStatsParentDir}
        biscuit cinread \
            {params.ref} {wildcards.BaseName}.bam {params.cinread_params}\
            -o {wildcards.BaseName}.biscuit.cinread.txt
        """
