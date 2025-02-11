configfile: "config/runtime_config.yaml"

rule astair_idbias:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        multiext("result/{fname}/{trimmer}/{aligner}/{deduper}/"
                 "{fname}{bqsr}",
                 "astair.IDbias_abundance_10bp_mod_site.pdf",
                 "astair.IDbias_abundance.pdf",
                 "astair.IDbias_indel_rate.pdf",
                 "astair.IDbias_mod_co-localize.pdf",
                 "astair.IDbias.stats")
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.fname.split('_')[1]],
        extra_params = (config["astair"]["idbias"]["extra_params"]
                        if config["astair"]["idbias"]["extra_params"] else ""),
        pattern      = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        astair idbias \
            -i {wildcards.fname}{params.pattern}.bam \
            -f {params.ref} \
            -d . \
            -l 150 -p -ni . \
            -t {threads} {params.extra_params}
        
        mv \
            {wildcards.fname}{params.pattern}_ID-bias_abundance_10bp_mod_site_plot.pdf \
            {wildcards.fname}{params.pattern}.astair.IDbias_abundance_10bp_mod_site.pdf
        mv \
            {wildcards.fname}{params.pattern}_ID-bias_abundance_plot.pdf \
            {wildcards.fname}{params.pattern}.astair.IDbias_abundance.pdf
        mv \
            {wildcards.fname}{params.pattern}_ID-bias_indel_rate_plot.pdf \
            {wildcards.fname}{params.pattern}.astair.IDbias_indel_rate.pdf
        mv \
            {wildcards.fname}{params.pattern}_ID-bias_modification_colocalisation_plot.pdf \
            {wildcards.fname}{params.pattern}.astair.IDbias_mod_co-localize.pdf
        mv \
            {wildcards.fname}{params.pattern}_ID-bias.stats \
            {wildcards.fname}{params.pattern}.astair.IDbias.stats
        """

rule astair_mbias:
    input:
        "result/{fname}/{trimmer}/{aligner}/{deduper}/{fname}{bqsr}bam"
    output:
        multiext("result/{fname}/{trimmer}/{aligner}/{deduper}/"
                 "{fname}{bqsr}",
                 "astair.Mbias.stats",
                 "astair.Mbias.pdf")
    params:
        ref          = lambda wildcards: config["ref"]["bwa-mem"][wildcards.fname.split('_')[1]],
        extra_params = (config["astair"]["mbias"]["extra_params"]
                        if config["astair"]["mbias"]["extra_params"] else ""),
        pattern      = lambda wildcards: "" if wildcards.bqsr == "." else ".bqsr"
    threads: 8
    conda:
        "conda.yaml"
    shell:
        """
        cd result/{wildcards.fname}/{wildcards.trimmer}/{wildcards.aligner}/{wildcards.deduper}
        astair mbias \
            -i {wildcards.fname}{params.pattern}.bam \
            -f {params.ref} -d . -l 150 -p -ni . \
            -t {threads} {params.extra_params}
        
        mv \
            {wildcards.fname}{params.pattern}_Mbias.txt \
            {wildcards.fname}{params.pattern}.astair.Mbias.stats
        mv \
            {wildcards.fname}{params.pattern}_M-bias_plot.pdf \
            {wildcards.fname}{params.pattern}.astair.Mbias.pdf
        """
