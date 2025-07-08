configfile: "config/runtime_config.yaml"
from textwrap import dedent

"""
epibed: 10 columns
chrom   start   end base    type    dinucleotide  penta-nucleotide  beta    m   u
chr1    790161  790162  G       CHH     CC      TTCCA   100     5       0

For TAPS, we need to flip the m and u columns.
"""

rule biscuit_count:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        vcf            = "result/{BaseName}/{CountParentDir}/biscuit/{BaseName}.vcf.gz",
        vcf_tbi        = "result/{BaseName}/{CountParentDir}/biscuit/{BaseName}.vcf.gz.tbi",
        epibed         = "result/{BaseName}/{CountParentDir}/biscuit/{BaseName}.epibed.gz"
    benchmark:
        "result/{BaseName}/{CountParentDir}/biscuit/{BaseName}.count.benchmark"
    params:
        ref            = lambda wildcards: config["ref"]["biscuit"][wildcards.BaseName.split('_')[1]],
        method         = lambda wildcards: wildcards.BaseName.split('_')[0][:2],
        pileup_params  = config["biscuit"]["pileup"]["extra_params"] or "",
        vcf2bed_params = config["biscuit"]["vcf2bed"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        export tmp_dir="result/{wildcards.BaseName}/{wildcards.CountParentDir}"
        mkdir -p ${{tmp_dir}}/biscuit
        biscuit pileup \\
          -@ {threads} {params.pileup_params} \\
          {params.ref} ${{tmp_dir}}/{wildcards.BaseName}.bam |\\
        bgzip -@ {threads} -o {output.vcf}
        tabix -p vcf {output.vcf}

        if [ {params.method} = "PS" ]; then
          biscuit vcf2bed \\
          {params.vcf2bed_params} {output.vcf} |\\
          awk 'BEGIN{{FS=OFS="\\t"}} {{
            if ($4 == "C") {{
                $4 = "+"
            }} else if ($4 == "G") {{
                $4 = "-"
            }} else {{
                $4 = "."
            }}
            print $1, $2, $3, $4, $5, $6, $7, $10, $9}}' |\\
          pigz --best > {output.epibed}
        else
          biscuit vcf2bed \\
          {params.vcf2bed_params} {output.vcf} |\\
          awk 'BEGIN{{FS=OFS="\\t"}} {{
            if ($4 == "C") {{
                $4 = "+"
            }} else if ($4 == "G") {{
                $4 = "-"
            }} else {{
                $4 = "."
            }}
            print $1, $2, $3, $4, $5, $6, $7, $9, $10}}' |\\
          pigz --best > {output.epibed}
        fi
        """)

