configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule haplotypecaller_call:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        gvcf                   = temp("result/{BaseName}/{CountParentDir}/haplotypecaller/{BaseName}.gvcf.gz"),
        vcf                    = "result/{BaseName}/{CountParentDir}/haplotypecaller/{BaseName}.vcf.gz"
    benchmark:
        "result/{BaseName}/{CountParentDir}/haplotypecaller/{BaseName}.hp_call.benchmark"
    params:
        ref                    = lambda wildcards: config["ref"]["gatk"]["fasta"][wildcards.BaseName.split('_')[1]],
        max_mem                = config["gatk"]["max_mem"],
        dbsnp_file             = config["gatk"]["haplotypecaller"]["dbsnp_file"],
        emit_ref_confidence    = config["gatk"]["haplotypecaller"]["emit_ref_confidence"],
        gam                    = config["gatk"]["haplotypecaller"]["genotype_assignment_method"],
        haplotypecaller_params = config["gatk"]["haplotypecaller"]["extra_params"] or "",
        genotypegvcfs_params   = config["gatk"]["genotypegvcfs"]["extra_params"] or ""
    threads: 64
    conda:
        "conda.yaml"
    shell:
        dedent("""
        export tmp_dir="result/{wildcards.BaseName}/{wildcards.CountParentDir}/haplotypecaller"
        gatk \\
          --java-options "-Xmx{params.max_mem}g -XX:ParallelGCThreads={threads}" HaplotypeCaller \\
          -R {params.ref} \\
          -I {input} \\
          -O {output.gvcf} \\
          --dbsnp {params.dbsnp_file} \\
          -gam {params.gam} \\
          --native-pair-hmm-threads {threads} \\
          -ERC {params.emit_ref_confidence} {params.haplotypecaller_params}
        gatk \\
          --java-options "-Xmx{params.max_mem}g -XX:ParallelGCThreads={threads}" GenotypeGVCFs \\
          -R {params.ref} \\
          -I {input} \\
          -V {output.gvcf} \\
          -O {output.vcf} \\
          --dbsnp {params.dbsnp_file} \\
          -gam {params.gam} {params.genotypegvcfs_params}
        """)


rule haplotypecaller_extract:
    input:
        "result/{BaseName}/{CountParentDir}/haplotypecaller/{BaseName}.vcf.gz"
    output:
        "result/{BaseName}/{CountParentDir}/haplotypecaller/{BaseName}.bedgraph.zst"
    benchmark:
        "result/{BaseName}/{CountParentDir}/haplotypecaller/{BaseName}.hp_extract.benchmark"
    params:
        method = lambda wildcards: wildcards.BaseName.split('_')[0][: 2],
    conda:
        "conda.yaml"
    threads: 8
    shell:
        dedent("""
        if [ "{params.method}" = "PS" ]; then
          m_col=2
          u_col=1
        else
          m_col=1
          u_col=2
        fi

        bcftools filter -i '
          ((REF == "C" && ALT == "T") || (REF == "G" && ALT == "A")) &&
          INFO/DP >= 5 &&
          INFO/MQ >= 10 &&
          INFO/QD > 2 &&
          ID == "."' {input} |\\
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%INFO/DP\\t[ %AD]\\n' |\\
        awk -v m_col=$m_col -v u_col=$u_col '
          BEGIN{{
          OFS="\\t"; FS="\\t";
          print "chrom", "start", "end", "strand", "depth", "m", "u";
          }}
          {{
            start= $2 - 1;

            if ($3 == "C") {{
              strand = "+"
            }} else {{ strand = "-" }}

            split($5, ad, ",");
            m = ad[m_col];
            u = ad[u_col];

            print $1, start, $2, strand, $4, m, u
          }}' | zstd -19 -T{threads} > {output}

        """)


rule methylationtypecaller_caller:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        "result/{BaseName}/{CountParentDir}/methylationtypecaller/{BaseName}.vcf.gz"
    benchmark:
        "result/{BaseName}/{CountParentDir}/methylationtypecaller/{BaseName}.mt_call.benchmark"
    params:
        ref      = lambda wildcards: config["ref"]["gatk"]["fasta"][wildcards.BaseName.split('_')[1]],
        ref_dict = lambda wildcards: config["ref"]["gatk"]["dict"][wildcards.BaseName.split('_')[1]],
        max_mem  = config["gatk"]["max_mem"]
    conda:
        "conda.yaml"
    threads: 8
    shell:
        dedent("""
        export OMP_NUM_THREADS={threads}
        export tmp_dir="result/{wildcards.BaseName}/{wildcards.CountParentDir}/methylationtypecaller"
        mkdir -p ${{tmp_dir}}
        gatk \\
          --java-options "-Xmx{params.max_mem}g -XX:ParallelGCThreads={threads}" \\
          MethylationTypeCaller \\
          -R {params.ref} \\
          --sequence-dictionary {params.ref_dict} \\
          -I {input} \\
          -O {output}
        """)


rule methylationtypecaller_extract:
    input:
        "result/{BaseName}/{CountParentDir}/methylationtypecaller/{BaseName}.vcf.gz"
    output:
        "result/{BaseName}/{CountParentDir}/methylationtypecaller/{BaseName}.bedgraph.zst"
    benchmark:
        "result/{BaseName}/{CountParentDir}/methylationtypecaller/{BaseName}.mt_extract.benchmark"
    params:
        method = lambda wildcards: wildcards.BaseName.split('_')[0][: 2]
    conda:
        "conda.yaml"
    threads: 8
    shell:
        dedent("""
        if [ "{params.method}" = "PS" ]; then
          m_col=5
          u_col=6
        else
          m_col=6
          u_col=5
        fi

        vcftools \\
          --gzvcf {input} \\
          --get-INFO CONVERTED_BASE_COV \\
          --get-INFO UNCONVERTED_BASE_COV \\
          --get-INFO REFERENCE_CONTEXT \\
          --stdout |\\
        awk -v m_col="$m_col" \\
            -v u_col="$u_col" 'BEGIN {{
          OFS="\\t";
          print "chrom", "start", "end", "strand", "context", "depth", "m", "u"
        }}
          NR == 1 {{next}}
          {{
            m = $m_col;
            u = $u_col;
            strand = ($3 == "C" ? "+" : "-");
            depth = m + u;
            print $1, $2 - 1, $2, strand, $7, depth, m, u
          }}' | zstd -19 -T{threads} > {output}
        """)
