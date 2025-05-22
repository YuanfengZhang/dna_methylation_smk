configfile: "config/runtime_config.yaml"
from textwrap import dedent


rule trie_dedup:
    input:
        r1            = "result/{BaseName}/{DedupParentDir}/{BaseName}.R1.fq.gz",
        r2            = "result/{BaseName}/{DedupParentDir}/{BaseName}.R2.fq.gz"
    output:
        tmp_inr1      = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.R1.fq"),
        tmp_inr2      = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/{BaseName}.R2.fq"),
        tmp_R1        = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/uniq_R1.txt"),
        tmp_R2        = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/uniq_R2.txt"),
        tmp_intersect = temp("result/{BaseName}/{DedupParentDir}/trie-dedup/uniq_intersect.txt"),
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
        pigz -dc {input.r1} > {output.tmp_inr1}
        python3 resources/TrieDedup-main/Python/TrieDedup.py -i {output.tmp_inr1} -v {params.extra_params} > {output.tmp_R1} 
        pigz -dc {input.r2} > {output.tmp_inr2}
        python3 resources/TrieDedup-main/Python/TrieDedup.py -i {output.tmp_inr2} -v {params.extra_params} > {output.tmp_R2} 
        comm -12 <(sort {output.tmp_R1}) <(sort {output.tmp_R2}) > {output.tmp_intersect}
        seqtk subseq {output.tmp_inr1} {output.tmp_intersect} | gzip -c > {output.r1}
        seqtk subseq {output.tmp_inr2} {output.tmp_intersect} | gzip -c > {output.r2}
        """)
