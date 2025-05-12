configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule rastair_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        "result/{BaseName}/{CountParentDir}/rastair/{BaseName}.rastair.mods.gz"
    benchmark:
        "result/{BaseName}/{CountParentDir}/rastair/{BaseName}.count.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["samtools"][wildcards.BaseName.split('_')[1]],
        mapq         = lambda wildcards: config["rastair"]["mapq"],
        base_qual    = lambda wildcards: config["rastair"]["base_quality"],
        error_model  = lambda wildcards: config["rastair"]["error_model"],
        extra_params = config["rastair"]["count"]["extra_params"] or ""
    threads: 8
    conda: "../astair/conda.yaml"  # use the pigz in that env.
    shell:
        dedent("""
        export PATH="$(pwd)/resources/rastair/target/x86_64-unknown-linux-gnu/release:$PATH"
        rastair call \\
            --fasta-file {params.ref} \\
            --min-mapq {params.mapq} \\
            --min-baseq {params.base_qual} \\
            --error-model {params.error_model} {params.extra_params} \\
            --threads {threads} --read-threads {threads} \\
            {input} | pigz -p {threads} --best > {output}
        """)
