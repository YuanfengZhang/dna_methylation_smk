configfile: "config/runtime_config.yaml"
from textwrap import dedent


# ! It is necessary to change the permission of output files after deepvariant,
# ! otherwise snakemake will not be able to collect the output files,
# ! which will lead to an fatal error.

# ! singularity 3.* does not support --cpus param. Delete line 84 if using singularity 3.*.

# no --threads param for bcftools filter and query because it only helps fast compression rather than calculating.


rule deepvariant_call:
    input:
        "result/{BaseName}/{CountParentDir}/{BaseName}.bam"
    output:
        vcf          = "result/{BaseName}/{CountParentDir}/deepvariant/{BaseName}.deepvariant.vcf.gz",
        gvcf         = "result/{BaseName}/{CountParentDir}/deepvariant/{BaseName}.deepvariant.gvcf.gz"
    benchmark:
        "result/{BaseName}/{CountParentDir}/deepvariant/{BaseName}.dv_call.benchmark"
    params:
        gpu_enabled  = "1" if config["deepvariant"]["gpu_enabled"] else "",
        container    = config["deepvariant"]["container"],
        docker_sig   = config["deepvariant"]["docker_sig"],
        ref_dir      = lambda wildcards: config["ref"]["deepvariant"]["ref_dir"][wildcards.BaseName.split('_')[1]],
        ref_fname    = lambda wildcards: config["ref"]["deepvariant"]["ref_name"][wildcards.BaseName.split('_')[1]],
        extra_params = config["deepvariant"]["extra_params"] or ""
    threads: 8
    shell:
        dedent("""
        export tmp_dir="$(pwd)/result/{wildcards.BaseName}/{wildcards.CountParentDir}"

        mkdir -p ${{tmp_dir}}/deepvariant/tmp
        
        # get correct dir containing the bam file
        input_bam="${{tmp_dir}}/{wildcards.BaseName}.bam"
        if [ -f "${{input_bam}}" ] && [ ! -L "${{input_bam}}" ]; then
          input_dir="${{tmp_dir}}"
        else
          parent_dir="$(dirname ${{tmp_dir}})"
          input_bam="${{parent_dir}}/{wildcards.BaseName}.bam"
          if [ -f "${{input_bam}}" ] && [ ! -L "${{input_bam}}" ]; then
            input_dir="${{parent_dir}}"
          else
            parent2_dir="$(dirname ${{parent_dir}})"
            input_bam="${{parent2_dir}}/{wildcards.BaseName}.bam"
            if [ -f "${{input_bam}}" ] && [ ! -L "${{input_bam}}" ]; then
              input_dir="${{parent2_dir}}"
            else
              echo "cannot find any directory containing {wildcards.BaseName}.bam!"
              echo "current_dir: ${{tmp_dir}}"
              exit 1
            fi
          fi
        fi

        if [ {params.docker_sig} == "docker" ]; then
          mount_flag="-v"
          if [ -z {params.gpu_enabled} ]; then
            container_run="docker run --rm"
          else
            container_run="docker run --rm --gpus 1"
          fi
        else
          mount_flag="-B"
          if [ -z {params.gpu_enabled} ]; then
            container_run="singularity exec"
          else
            container_run="singularity exec --nv"
          fi
        fi

        echo -e "container run method:\\t${{container_run}}"
        echo -e "container:\\t{params.container}"
        echo -e "input dir:\\t${{input_dir}}"
        echo -e "output dir:\\t${{tmp_dir}}/deepvariant"
        echo -e "ref dir:\\t{params.ref_dir}"
        echo -e "mount flag:\\t${{mount_flag}}"

        ${{container_run}} \\
          ${{mount_flag}} {params.ref_dir}:/data/ref \\
          ${{mount_flag}} $(pwd)/resources:/data/resources \\
          ${{mount_flag}} ${{input_dir}}:/data/input \\
          ${{mount_flag}} ${{tmp_dir}}/deepvariant:/data/output \\
          ${{mount_flag}} ${{tmp_dir}}/deepvariant/tmp:/tmp \\
          --cpus {threads} \\
          {params.container} \\
          /opt/deepvariant/bin/run_deepvariant \\
          --model_type=WGS \\
          --ref=/data/ref/{params.ref_fname} \\
          --reads=/data/input/{wildcards.BaseName}.bam \\
          --output_vcf=/data/output/{wildcards.BaseName}.deepvariant.vcf.gz \\
          --output_gvcf=/data/output/{wildcards.BaseName}.deepvariant.gvcf.gz \\
          --num_shards={threads} \\
          --vcf_stats_report=true \\
          --logging_dir=/data/output/logs \\
          --haploid_contigs="chrX,chrY" \\
          --par_regions_bed="/data/resources/GRCh38_PAR.bed" \\
          --dry_run=false && chmod 777 -R /data/output
          """)


rule deepvariant_extract:
    input:
        "result/{BaseName}/{CountParentDir}/deepvariant/{BaseName}.deepvariant.vcf.gz"
    output:
        "result/{BaseName}/{CountParentDir}/deepvariant/{BaseName}.deepvariant.bedgraph.zst"
    benchmark:
        "result/{BaseName}/{CountParentDir}/deepvariant/{BaseName}.dv_extract.benchmark"
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

        bcftools filter \\
          -i '((REF=="C" & ALT=="T") | (REF=="G" & ALT=="A")) & (FILTER=="RefCall" | FILTER=="PASS")' \\
          {input} |\\
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t[ %DP]\\t[ %AD]\\n' |\\
        awk -v m_col="$m_col" -v u_col="$u_col" 'BEGIN {{
          OFS = "\\t";
          print "chrom", "start", "end", "strand", "depth", "m", "u"
        }}
        {{
          split($5, ad, ",")
          m = ad[m_col]
          u = ad[u_col]
          strand = ($3 == "C" ? "+" : "-");
          print $1, $2 - 1, $2, strand, $4, m, u
        }}' | zstd -19 -T{threads} > {output}
        """)
