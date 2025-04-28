configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule astair_idbias:
    input:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.bam"
    output:
        multiext("result/{BaseName}/{BamStatsParentDir}/{BaseName}",
                 ".astair.IDbias_abundance_10bp_mod_site.pdf",
                 ".astair.IDbias_abundance.pdf",
                 ".astair.IDbias_indel_rate.pdf",
                 ".astair.IDbias_mod_co-localize.pdf",
                 ".astair.IDbias.stats")
    benchmark:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.idbias.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        extra_params = config["astair"]["idbias"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.BamStatsParentDir}
        astair idbias \\
            -i {wildcards.BaseName}.bam \\
            -f {params.ref} \\
            -d . \\
            -l 150 -p -ni . \\
            -t {threads} {params.extra_params}
        
        mv \\
            {wildcards.BaseName}_ID-bias_abundance_10bp_mod_site_plot.pdf \\
            {wildcards.BaseName}.astair.IDbias_abundance_10bp_mod_site.pdf
        mv \\
            {wildcards.BaseName}_ID-bias_abundance_plot.pdf \\
            {wildcards.BaseName}.astair.IDbias_abundance.pdf
        mv \\
            {wildcards.BaseName}_ID-bias_indel_rate_plot.pdf \\
            {wildcards.BaseName}.astair.IDbias_indel_rate.pdf
        mv \\
            {wildcards.BaseName}_ID-bias_modification_colocalisation_plot.pdf \\
            {wildcards.BaseName}.astair.IDbias_mod_co-localize.pdf
        mv \\
            {wildcards.BaseName}_ID-bias.stats \\
            {wildcards.BaseName}.astair.IDbias.stats
        """)

rule astair_mbias:
    input:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.bam"
    output:
        multiext("result/{BaseName}/{BamStatsParentDir}/"
                 "{BaseName}",
                 ".astair.Mbias.stats",
                 ".astair.Mbias.pdf")
    benchmark:
        "result/{BaseName}/{BamStatsParentDir}/{BaseName}.mbias.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.BaseName.split('_')[1]],
        extra_params = config["astair"]["mbias"]["extra_params"]
                        if config["astair"]["mbias"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        cd result/{wildcards.BaseName}/{wildcards.BamStatsParentDir}
        astair mbias \\
            -i {wildcards.BaseName}.bam \\
            -f {params.ref} -d . -l 150 -p -ni . \\
            -t {threads} {params.extra_params}
        
        mv \\
            {wildcards.BaseName}_Mbias.txt \\
            {wildcards.BaseName}.astair.Mbias.stats
        mv \\
            {wildcards.BaseName}_M-bias_plot.pdf \\
            {wildcards.BaseName}.astair.Mbias.pdf
        """)
