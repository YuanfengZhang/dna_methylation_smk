"""
mkdir -p ./BL/methylGrapher
vg construct -r BL.fa -t 64 -p > BL.vg
vg view BL.vg > BL.gfa
(bioconda::vg)
methylGrapher PrepareGenome -gfa ./BL.gfa -prefix ./BL -t 24
methylGrapher Main \
    -t 24 \
    -work_dir /mnt/eqa/zhangyuanfeng/methylation/dna_methylation_smk/test/tmp \
    -index_prefix /hot_warm_data/ref/quartet/DNA/custom_genome/BL/methylGrapher/BL \
    -fq1 /mnt/eqa/zhangyuanfeng/methylation/dna_methylation_smk/input/BS_HF_1_MX.R1.fq.gz \
    -fq2 /mnt/eqa/zhangyuanfeng/methylation/dna_methylation_smk/input/BS_HF_1_MX.R2.fq.gz
"""