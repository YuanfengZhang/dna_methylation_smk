configfile: "config/runtime_config.yaml"

"""
rule bismark_extractor will generate these files:
├── CHG_OB_${fname}.txt.gz -> ${fname}.CHG_OB.txt.gz
├── CHG_OT_${fname}.txt.gz -> ${fname}.CHG_OT.txt.gz
├── CHH_OB_${fname}.txt.gz -> ${fname}.CHH_OB.txt.gz
├── CHH_OT_${fname}.txt.gz -> ${fname}.CHH_OT.txt.gz
├── CpG_OB_${fname}.txt.gz -> ${fname}.CpG_OB.txt.gz
├── CpG_OT_${fname}.txt.gz -> ${fname}.CpG_OT.txt.gz
├── ${fname}.bedGraph.gz -> ${fname}.bedgraph.gz (renamed)
├── ${fname}.bedGraph.gz.bismark.zero.cov
├── ${fname}.bedGraph_UCSC.bedGraph.gz -> ${fname}.ucsc.bedgraph.gz (renamed)
├── ${fname}.bismark.cov.gz
├── ${fname}.M-bias_R1.png
├── ${fname}.M-bias_R2.png
├── ${fname}.M-bias.txt
└── ${fname}.splitting_report.txt

or

├── CHG_OB_${fname}.bqsr.txt.gz -> ${fname}.bqsr.CHG_OB.txt.gz
├── CHG_OT_${fname}.bqsr.txt.gz -> ${fname}.bqsr.CHG_OT.txt.gz
├── CHH_OB_${fname}.bqsr.txt.gz -> ${fname}.bqsr.CHH_OT.txt.gz
├── CHH_OT_${fname}.bqsr.txt.gz -> ${fname}.bqsr.CHH_OT.txt.gz
├── CpG_OB_${fname}.bqsr.txt.gz -> ${fname}.bqsr.CpG_OB.txt.gz
├── CpG_OT_${fname}.bqsr.txt.gz -> ${fname}.bqsr.CpG_OT.txt.gz
├── ${fname}.bqsr.bedGraph.gz -> ${fname}.bqsr.bedgraph.gz (renamed)
├── ${fname}.bqsr.bedGraph.gz.bismark.zero.cov
├── ${fname}.bqsr.bedGraph_UCSC.bedGraph.gz -> ${fname}.bqsr.ucsc.bedgraph.gz (renamed)
├── ${fname}.bqsr.bismark.cov.gz
├── ${fname}.bqsr.M-bias_R1.png
├── ${fname}.bqsr.M-bias_R2.png
├── ${fname}.bqsr.M-bias.txt
└── ${fname}.bqsr.splitting_report.txt

depends on whether the RECALIBRATE is True or False.
"""
rule bismark_extractor:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        multiext("result/{fname}/{trimmer}/{aligner}/{deduper}"
                 "/bismark/{fname}{bqsr}",
                 "bedgraph.gz", "ucsc_bedgraph.gz",
                 "CHG_OB.txt.gz", "CHG_OT.txt.gz",
                 "CHH_OB.txt.gz", "CHH_OT.txt.gz",
                 "CpG_OB.txt.gz", "CpG_OT.txt.gz",
                 "bismark.cov.gz", "splitting_report",
                 "bismark.zero.cov", "M-bias.txt",
                 "M-bias_R1.png", "M-bias_R2.png")
    params:
        ref           = lambda wildcards: config["ref"]["bismark-bowtie2"][wildcards.fname.split('_')[1]],
        extra_params  = (config["bismark"]["exactor"]["extra_params"]
                         if config["bismark"]["exactor"]["extra_params"] else ""),
        pattern       = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        samtools sort \
            -n -@ {threads} \
            -o {wildcards.fname}{params.pattern}.sortn.bam \
            {wildcards.fname}{params.pattern}.bam
        bismark_methylation_extractor \
            {wildcards.fname}{params.pattern}.sortn.bam \
            -o bismark \
            --genome_folder {params.ref} \
            --gzip --bedGraph --cutoff 5 --zero_based --scaffolds \
            --parallel {threads} --buffer_size 20G --ucsc \
            --paired-end --no_header \
            {params.extra_params}
        cd bismark
        mv \
            {wildcards.fname}{params.pattern}.sortn.bedGraph.gz \
            {wildcards.fname}{params.pattern}.bedgraph.gz
        mv \
            {wildcards.fname}{params.pattern}.sortn.bedGraph_UCSC.bedGraph.gz \
            {wildcards.fname}{params.pattern}.ucsc_bedgraph.gz

        mv \
            CHG_OB_{wildcards.fname}{params.pattern}.sortn.txt.gz \
            {wildcards.fname}{params.pattern}.CHG_OB.txt.gz
        mv \
            CHG_OT_{wildcards.fname}{params.pattern}.sortn.txt.gz \
            {wildcards.fname}{params.pattern}.CHG_OT.txt.gz
        mv \
            CHH_OB_{wildcards.fname}{params.pattern}.sortn.txt.gz \
            {wildcards.fname}{params.pattern}.CHH_OB.txt.gz
        mv \
            CHH_OT_{wildcards.fname}{params.pattern}.sortn.txt.gz \
            {wildcards.fname}{params.pattern}.CHH_OT.txt.gz
        mv \
            CpG_OB_{wildcards.fname}{params.pattern}.sortn.txt.gz \
            {wildcards.fname}{params.pattern}.CpG_OB.txt.gz
        mv \
            CpG_OT_{wildcards.fname}{params.pattern}.sortn.txt.gz \
            {wildcards.fname}{params.pattern}.CpG_OT.txt.gz
        mv \
            {wildcards.fname}{params.pattern}.sortn_splitting_report.txt \
            {wildcards.fname}{params.pattern}.splitting_report
        mv \
            {wildcards.fname}{params.pattern}.sortn.bismark.cov.gz \
            {wildcards.fname}{params.pattern}.bismark.cov.gz
        mv \
            {wildcards.fname}{params.pattern}.sortn.bedGraph.gz.bismark.zero.cov \
            {wildcards.fname}{params.pattern}.bismark.zero.cov
        mv \
            {wildcards.fname}{params.pattern}.sortn.M-bias_R1.png \
            {wildcards.fname}{params.pattern}.M-bias_R1.png
        mv \
            {wildcards.fname}{params.pattern}.sortn.M-bias_R2.png \
            {wildcards.fname}{params.pattern}.M-bias_R2.png
        mv \
            {wildcards.fname}{params.pattern}.sortn.M-bias.txt \
            {wildcards.fname}{params.pattern}.M-bias.txt
        """

"""
rule bismark_c2c will generate these files:

├── ${fname}.CpG_report.merged_CpG_evidence.cov.gz -> ${fname}.c2c.cov.gz (renamed)
├── ${fname}.CpG_report.txt.gz -> ${fname}.c2c.report.gz (renamed)
└── ${fname}.cytosine_context_summary.txt -> ${fname}.c2c.summary (renamed)

or

├── ${fname}.bqsr.CpG_report.merged_CpG_evidence.cov.gz -> ${fname}.bqsr.c2c.cov.gz (renamed)
├── ${fname}.bqsr.CpG_report.txt.gz -> ${fname}.bqsr.c2c.report.gz (renamed)
└── ${fname}.bqsr.cytosine_context_summary.txt -> ${fname}.bqsr.c2c.summary (renamed)
"""
rule bismark_c2c:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        multiext("result/{fname}/{trimmer}/{aligner}/{deduper}"
                 "/bismark/{fname}{bqsr}",
                 "c2c.cov.gz",
                 "c2c.report.gz",
                 "c2c.summary")
    params:
        ref           = lambda wildcards: config["ref"]["bismark-bowtie2"][wildcards.fname.split('_')[1]],
        extra_params  = (config["bismark"]["coverage2cytosine"]["extra_params"]
                         if config["bismark"]["coverage2cytosine"]["extra_params"] else ""),
        pattern       = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        coverage2cytosine \
            {wildcards.fname}{params.pattern}.bam \
            --output bismark/{wildcards.fname}{params.pattern} \
            --genome_folder {params.ref} \
            --merge_CpG --ff --zero_based --gzip \
            {params.extra_params}
        cd bismark
        mv \
            {wildcards.fname}{params.pattern}.CpG_report.merged_CpG_evidence.cov.gz \
            {wildcards.fname}{params.pattern}.c2c.cov.gz
        mv \
            {wildcards.fname}{params.pattern}.CpG_report.txt.gz \
            {wildcards.fname}{params.pattern}.c2c.report.gz
        mv \
            {wildcards.fname}{params.pattern}.cytosine_context_summary.txt \
            {wildcards.fname}{params.pattern}.c2c.summary
        """
