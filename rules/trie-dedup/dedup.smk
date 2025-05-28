configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule trie_dedup:
    input:
        r1            = "result/{BaseName}/{DedupParentDir}/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{DedupParentDir}/{BaseName}.R2.fq.gz"
    output:
        uniq_read_id  = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.unique_read_list")
        merged_fq     = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.merged.fq"),
        subseq_fq     = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.subseq.fq"),
        r1            = "result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.R2.fq.gz"
    benchmark:
        "result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.dedup.benchmark"
    params:
        method        = config["trie-dedup"]["dedup"]["method"], # sortuniq / trie / pairwise
        extra_params  = config["trie-dedup"]["dedup"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        mkdir -p "result/{wildcards.BaseName}/{wildcards.DedupParentDir}/trie-dedup"
    
        seqtk mergepe {input.r1} {input.r2} |\\
          seqtk seq -q 10 -n N - > {output.merged_fq}
        
        resources/TrieDedup/Cpp/bin/TrieDedup \\
          {params.method} -f readID \\
          -o {output.uniq_read_id} {output.merged_fq}

        seqtk subseq {output.tmp_fq} {output.uniq_read_id} > {output.subseq_fq}
        grep '^@.*/1$' -A 3 --no-group-separator | pigz --best -p {threads} > {output.r1}
        grep '^@.*/2$' -A 3 --no-group-separator | pigz --best -p {threads} > {output.r2}
        """)
