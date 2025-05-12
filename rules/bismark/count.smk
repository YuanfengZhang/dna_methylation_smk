configfile: "config/runtime_config.yaml"
from textwrap import dedent


"""
rule bismark_extractor will generate these files:
├── CHG_OB_${BaseName}.txt.gz -> ${BaseName}.CHG_OB.txt.gz
├── CHG_OT_${BaseName}.txt.gz -> ${BaseName}.CHG_OT.txt.gz
├── CHH_OB_${BaseName}.txt.gz -> ${BaseName}.CHH_OB.txt.gz
├── CHH_OT_${BaseName}.txt.gz -> ${BaseName}.CHH_OT.txt.gz
├── CpG_OB_${BaseName}.txt.gz -> ${BaseName}.CpG_OB.txt.gz
├── CpG_OT_${BaseName}.txt.gz -> ${BaseName}.CpG_OT.txt.gz
├── ${BaseName}.bedGraph.gz -> ${BaseName}.bedgraph.gz (renamed)
├── ${BaseName}.bedGraph.gz.bismark.zero.cov
├── ${BaseName}.bedGraph_UCSC.bedGraph.gz -> ${BaseName}.ucsc.bedgraph.gz (renamed)
├── ${BaseName}.bismark.cov.gz
├── ${BaseName}.M-bias_R1.png
├── ${BaseName}.M-bias_R2.png
├── ${BaseName}.M-bias.txt
└── ${BaseName}.splitting_report.txt
"""

rule bismark_extractor:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        multiext("result/{BaseName}/{CountParentDir}"
                 "/bismark/{BaseName}.",
                 "bedgraph.gz", "ucsc_bedgraph.gz",
                 "CHG_OB.txt.gz", "CHG_OT.txt.gz",
                 "CHH_OB.txt.gz", "CHH_OT.txt.gz",
                 "CpG_OB.txt.gz", "CpG_OT.txt.gz",
                 "bismark.cov.gz", "splitting_report",
                 "bismark.zero.cov", "M-bias.txt",
                 "M-bias_R1.png", "M-bias_R2.png")
    benchmark:
        "result/{BaseName}/{CountParentDir}/bismark/{BaseName}.extractor.benchmark"
    params:
        ref           = lambda wildcards: config["ref"]["bismark-bowtie2"][wildcards.BaseName.split('_')[1]],
        extra_params  = config["bismark"]["exactor"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.CountParentDir}
        samtools sort \\
            -n -@ {threads} \\
            -o {wildcards.BaseName}.sortn.bam \\
            {wildcards.BaseName}.bam
        bismark_methylation_extractor \\
            {wildcards.BaseName}.sortn.bam \\
            -o bismark \\
            --genome_folder {params.ref} \\
            --gzip --bedGraph --cutoff 5 --zero_based --scaffolds \\
            --parallel {threads} --buffer_size 20G --ucsc \\
            --paired-end --no_header \\
            {params.extra_params}
        cd bismark
        mv \\
            {wildcards.BaseName}.sortn.bedGraph.gz \\
            {wildcards.BaseName}.bedgraph.gz
        mv \\
            {wildcards.BaseName}.sortn.bedGraph_UCSC.bedGraph.gz \\
            {wildcards.BaseName}.ucsc_bedgraph.gz

        mv \\
            CHG_OB_{wildcards.BaseName}.sortn.txt.gz \\
            {wildcards.BaseName}.CHG_OB.txt.gz
        mv \\
            CHG_OT_{wildcards.BaseName}.sortn.txt.gz \\
            {wildcards.BaseName}.CHG_OT.txt.gz
        mv \\
            CHH_OB_{wildcards.BaseName}.sortn.txt.gz \\
            {wildcards.BaseName}.CHH_OB.txt.gz
        mv \\
            CHH_OT_{wildcards.BaseName}.sortn.txt.gz \\
            {wildcards.BaseName}.CHH_OT.txt.gz
        mv \\
            CpG_OB_{wildcards.BaseName}.sortn.txt.gz \\
            {wildcards.BaseName}.CpG_OB.txt.gz
        mv \\
            CpG_OT_{wildcards.BaseName}.sortn.txt.gz \\
            {wildcards.BaseName}.CpG_OT.txt.gz
        mv \\
            {wildcards.BaseName}.sortn_splitting_report.txt \\
            {wildcards.BaseName}.splitting_report
        mv \\
            {wildcards.BaseName}.sortn.bismark.cov.gz \\
            {wildcards.BaseName}.bismark.cov.gz
        mv \\
            {wildcards.BaseName}.sortn.bedGraph.gz.bismark.zero.cov \\
            {wildcards.BaseName}.bismark.zero.cov
        mv \\
            {wildcards.BaseName}.sortn.M-bias_R1.png \\
            {wildcards.BaseName}.M-bias_R1.png
        mv \\
            {wildcards.BaseName}.sortn.M-bias_R2.png \\
            {wildcards.BaseName}.M-bias_R2.png
        mv \\
            {wildcards.BaseName}.sortn.M-bias.txt \\
            {wildcards.BaseName}.M-bias.txt
        """)

"""
rule bismark_c2c will generate these files:

├── ${BaseName}.CpG_report.merged_CpG_evidence.cov.gz -> ${BaseName}.c2c.cov.gz (renamed)
├── ${BaseName}.CpG_report.txt.gz -> ${BaseName}.c2c.report.gz (renamed)
└── ${BaseName}.cytosine_context_summary.txt -> ${BaseName}.c2c.summary (renamed)
"""

rule bismark_c2c:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}/bismark/{BaseName}.bismark.cov.gz"
    output:
        multiext("result/{BaseName}/{CountParentDir}"
                 "/bismark/{BaseName}",
                 ".c2c.cov.gz",
                 ".c2c.report.gz",
                 ".c2c.summary")
    benchmark:
        "result/{BaseName}/{CountParentDir}/bismark/{BaseName}.c2c.benchmark"
    params:
        ref           = lambda wildcards: config["ref"]["bismark-bowtie2"][wildcards.BaseName.split('_')[1]],
        extra_params  = config["bismark"]["coverage2cytosine"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.CountParentDir}/bismark
        coverage2cytosine \\
            {wildcards.BaseName}.bismark.cov.gz \\
            --output {wildcards.BaseName} \\
            --genome_folder {params.ref} \\
            --merge_CpG --ff --zero_based --gzip \\
            {params.extra_params}
        cd bismark
        mv \\
            {wildcards.BaseName}.CpG_report.merged_CpG_evidence.cov.gz \\
            {wildcards.BaseName}.c2c.cov.gz
        mv \\
            {wildcards.BaseName}.CpG_report.txt.gz \\
            {wildcards.BaseName}.c2c.report.gz
        mv \\
            {wildcards.BaseName}.cytosine_context_summary.txt \\
            {wildcards.BaseName}.c2c.summary
        """)
