configfile: "config/runtime_config.yaml"
from textwrap import dedent


# to build the index:
# resources/FAME/FAME --genome ref.fa --store_index ref.fa.fame


# ! FAME reports the m+, m-, u+, u- counts in the output tsv file.
# ! We need to retrieve the sequence from fasta file and separate the strands.


rule fame_fq2bedgraph:
    input:
        r1           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R1.fq.gz",
        r2           = "result/{BaseName}/{AlignParentDir}/{BaseName}.R2.fq.gz"
    output:
        tsv          = temp("result/{BaseName}/{AlignParentDir}/fame/{BaseName}_cpg.tsv"),
        fame         = temp("result/{BaseName}/{AlignParentDir}/fame/{BaseName}.fame"),
        bedgraph     = "result/{BaseName}/{AlignParentDir}/fame/{BaseName}.bedgraph.zst"
    benchmark:
        "result/{BaseName}/{AlignParentDir}/fame/{BaseName}.count.benchmark"
    params:
        ref          = lambda wildcards: config["ref"]["fame"][wildcards.BaseName.split('_')[1]],
        extra_params = config["fame"]["extra_params"] or ""
    threads: 8
    conda: "conda.yaml"
    shell:
        dedent("""
        tmp_dir="result/{wildcards.BaseName}/{wildcards.AlignParentDir}/fame"
        mkdir -p ${{tmp_dir}}/tmp
        resources/FAME/FAME \\
          -r1 {input.r1} -r2 {input.r2} \\
          --load_index {params.ref}.fame {params.extra_params} \\
          -o ${{tmp_dir}}/{wildcards.BaseName}
        awk -F '\\t' '{{sum=$3+$4+$5+$6}} sum >= 5' {output.tsv} |\
        sort -k1,2n > {output.fame}


        awk -F'\\t' -v tmp_dir="$tmp_dir/tmp" '
          BEGIN {{
            OFS = "\\t";
            count = 0;
            out_file_count = 0;
        }}
        {{
            end = $2 + 2;

            if (count % 1000000 == 0) {{
                if (out_file_count > 0) {{
                    close(out_file);
                }}
                out_file_count++;
                out_file = tmp_dir "/" out_file_count ".bed";
            }}

            print $1, $2, end, $3, $4, $5, $6 > out_file;
            count++;
        }}' {output.fame}

        # for f in ${{tmp_dir}}/tmp/*bed; do
        #     fname=$(basename $f)
        #     bedtools getfasta \\
        #       -fi {params.ref} \\
        #       -bed $f \\
        #       -bedOut >$f.seq &
        # done
        # wait
        find ${{tmp_dir}}/tmp -type f -name "*.bed" | parallel -j {threads} '
          bedtools getfasta \\
            -fi {params.ref} \\
            -bed {{}} \\
            -bedOut > {{.}}.seq'

        cat ${{tmp_dir}}/tmp/*.seq | \\
        awk -F'\\t' '
        BEGIN {{
          print "chrom", "start", "end", "m", "u", "strand", "depth"
        }}
        {{
          if ($8 ~ /[Cc][Gg]/) {{
            depth1 = $4 + $5
            if (depth1 >= 5)
              print $1, $2, $2 + 1, $4, $5, "+", depth1

            depth2 = $6 + $7
            if (depth2 >= 5)
              print $1, $2 + 1, $2 + 2, $6, $7, "-", depth2
        }}
        else if ($8 ~ /[Gg][Cc]/) {{
          depth1 = $6 + $7
          if (depth1 >= 5)
              print $1, $2, $2 + 1, $6, $7, "+", depth1

          depth2 = $4 + $5
          if (depth2 >= 5)
              print $1, $2 + 1, $2 + 2, $4, $5, "-", depth2
        }}
        else {{
          next
        }}
    }}' OFS='\\t' |\\
        sort -k1,1 -k2,2n |\\
        zstd -19 -T{threads} > {output.bedgraph}
        """)
