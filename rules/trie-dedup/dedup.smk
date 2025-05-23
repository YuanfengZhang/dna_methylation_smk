configfile: "config/runtime_config.yaml"
from textwrap import dedent

# ! Do not use this.
# ! Since it does not support paired-end fastq files,
# ! the output files are troblesome.


rule trie_dedup:
    input:
        r1            = "result/{BaseName}/{DedupParentDir}/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{DedupParentDir}/{BaseName}.R2.fq.gz"
    output:
        tmp_inr1      = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.R1.fq"),
        tmp_inr2      = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.R2.fq"),
        tmp_id1       = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.R1.id"),
        tmp_id2       = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.R2.id"),
        r1            = "result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.R2.fq.gz",
    benchmark:
        "result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.dedup.benchmark"
    params:
        extra_params  = lambda wildcards: config["trie-dedup"]["dedup"]["extra_params"] or ""
    threads: 8
    conda:
        "conda.yaml"
    shell:
        dedent("""
        pigz -dc {input.r1} > {output.tmp_inr1} &
        pigz -dc {input.r2} > {output.tmp_inr2} &
        wait
        resources/TrieDedup/Cpp/bin/TrieDedup \
          -f readID -o {output.tmp_id1} {output.tmp_inr1}
        resources/TrieDedup/Cpp/bin/TrieDedup \
          -f readID -o {output.tmp_id2} {output.tmp_inr2}
        
        grep -f {output.tmp_id1} {output.tmp_inr1} | pigz -p {threads} > {output.r1}
        grep -f {output.tmp_id2} {output.tmp_inr2} | pigz -p {threads} > {output.r2}
        """)
