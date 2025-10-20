pairs=("D6_vs_D5" "D6_vs_F7" "D6_vs_M8" "BC_vs_BL")
for pair in "${pairs[@]}"; do
  python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_ref_prepare.py \
    pseudo_lab\
    -i /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/final/full_seq_info \
    -o /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc/pseudo_lab
done

pairs=("D6_vs_D5" "D6_vs_F7" "D6_vs_M8" "BC_vs_BL")
p_types=("p" "q")
dmc_root_dir="/mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc"
cytosine_truset_dir="/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/final/full_seq_info"
dmc_truset_dir="/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc"

for pair in "${pairs[@]}"; do
  for p_type in "${p_types[@]}"; do
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_ref_prepare.py \
      merge_1040 \
      -m methylkit --p-type ${p_type} \
      --bs0 "${dmc_root_dir}/methylkit/output/BS0_${pair}.methylkit" \
      --em0 "${dmc_root_dir}/methylkit/output/EM0_${pair}.methylkit" \
      -r "${cytosine_truset_dir}" \
      -o "${dmc_truset_dir}"
  done
done


# ! Before
labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
pairs=("D6,D5" "D6,F7" "D6,M8" "BC,BL")
tool="methylkit"
for lab in "${labs[@]}"; do
  for pair in "${pairs[@]}"; do
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_prepare.py \
      -i /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/merged \
      -o /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/before/c_dmc/${tool}/input \
      -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/final/full_seq_info \
      -l ${lab} -s ${pair} -m ${tool} --beta-col predicted_beta \
      -d 5  -c chr17 -v debug -a prepare
  done
done

labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
pairs=("D6,D5" "D6,F7" "D6,M8" "BC,BL")
tool="methylsig"
for lab in "${labs[@]}"; do
  for pair in "${pairs[@]}"; do
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_prepare.py \
      -i /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/merged \
      -o /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/before/c_dmc/${tool}/input \
      -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/final/full_seq_info \
      -l ${lab} -s ${pair} -m ${tool} --beta-col predicted_beta \
      -d 5  -c chr17 -v debug -a prepare
  done
done

labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
pairs=("D6,D5" "D6,F7" "D6,M8" "BC,BL")
tool="cpgtools"
for lab in "${labs[@]}"; do
  for pair in "${pairs[@]}"; do
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_prepare.py \
      -i /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/formatted \
      -o /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/before/c_dmc/${tool}/input \
      -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/final/full_seq_info \
      -l ${lab} -s ${pair} -m ${tool} --beta-col predicted_beta \
      -d 5  -c chr17 -v debug -a prepare
  done
done

labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
cpgtools_pyscript="/home/zhangyuanfeng/mambaforge/envs/cpgtools/bin/dmc_fisher.py"
cpgtools_dir="/mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/before/c_dmc/cpgtools"
for lab in "${labs[@]}"; do
  time python "${cpgtools_pyscript}" \
    -i ${cpgtools_dir}/input/${lab}_D6_vs_D5.bed \
    -g ${cpgtools_dir}/input/${lab}_D6_vs_D5.group \
    -o ${cpgtools_dir}/output/${lab}_D6_vs_D5.cpgtools &
  time python "${cpgtools_pyscript}" \
    -i ${cpgtools_dir}/input/${lab}_D6_vs_F7.bed \
    -g ${cpgtools_dir}/input/${lab}_D6_vs_F7.group \
    -o ${cpgtools_dir}/output/${lab}_D6_vs_F7.cpgtools &
  time python "${cpgtools_pyscript}" \
    -i ${cpgtools_dir}/input/${lab}_D6_vs_M8.bed \
    -g ${cpgtools_dir}/input/${lab}_D6_vs_M8.group \
    -o ${cpgtools_dir}/output/${lab}_D6_vs_M8.cpgtools &
  time python "${cpgtools_pyscript}" \
    -i ${cpgtools_dir}/input/${lab}_BC_vs_BL.bed \
    -g ${cpgtools_dir}/input/${lab}_BC_vs_BL.group \
    -o ${cpgtools_dir}/output/${lab}_BC_vs_BL.cpgtools & wait
done

labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
for lab in ${labs[@]}; do
  docker run --rm \
    -v /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/before/c_dmc/methylkit/input:/data/input \
    -v /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/before/c_dmc/methylkit/output:/data/output \
    -v /mnt/eqa/zhangyuanfeng/methylation/src/evaluation:/data/bin \
    methylkit:1.33.3 \
    Rscript /data/bin/methylkit.R -l ${lab}
done

labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
for lab in ${labs[@]}; do
  docker run --rm \
    -v /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/before/c_dmc/methylsig/input:/data/input \
    -v /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/before/c_dmc/methylsig/output:/data/output \
    -v /mnt/eqa/zhangyuanfeng/methylation/src/evaluation:/data/bin \
    methylsig:1.19.0 \
    Rscript /data/bin/methylsig.R -l ${lab}
done


# ! After
labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
pairs=("D6,D5" "D6,F7" "D6,M8" "BC,BL")
tool="methylkit"
for lab in "${labs[@]}"; do
  for pair in "${pairs[@]}"; do
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_prepare.py \
      -i /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/for_dmc \
      -o /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/after/c_dmc/${tool}/input \
      -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/final/full_seq_info \
      -l ${lab} -s ${pair} -m ${tool} --beta-col actual_beta \
      -d 5  -c chr17 -v debug -a prepare
  done
done

labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
pairs=("D6,D5" "D6,F7" "D6,M8" "BC,BL")
tool="methylsig"
for lab in "${labs[@]}"; do
  for pair in "${pairs[@]}"; do
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_prepare.py \
      -i /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/for_dmc \
      -o /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/after/c_dmc/${tool}/input \
      -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/final/full_seq_info \
      -l ${lab} -s ${pair} -m ${tool} --beta-col actual_beta \
      -d 5  -c chr17 -v debug -a prepare
  done
done

labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
pairs=("D6,D5" "D6,F7" "D6,M8" "BC,BL")
tool="cpgtools"
for lab in "${labs[@]}"; do
  for pair in "${pairs[@]}"; do
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_prepare.py \
      -i /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/for_dmc \
      -o /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/after/c_dmc/${tool}/input \
      -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/final/full_seq_info \
      -l ${lab} -s ${pair} -m ${tool} --beta-col actual_beta \
      -d 5  -c chr17 -v debug -a prepare
  done
done

labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
cpgtools_pyscript="/home/zhangyuanfeng/mambaforge/envs/cpgtools/bin/dmc_fisher.py"
cpgtools_dir="/mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/after/c_dmc/cpgtools"
for lab in "${labs[@]}"; do
  time python "${cpgtools_pyscript}" \
    -i ${cpgtools_dir}/input/${lab}_D6_vs_D5.bed \
    -g ${cpgtools_dir}/input/${lab}_D6_vs_D5.group \
    -o ${cpgtools_dir}/output/${lab}_D6_vs_D5.cpgtools &
  time python "${cpgtools_pyscript}" \
    -i ${cpgtools_dir}/input/${lab}_D6_vs_F7.bed \
    -g ${cpgtools_dir}/input/${lab}_D6_vs_F7.group \
    -o ${cpgtools_dir}/output/${lab}_D6_vs_F7.cpgtools &
  time python "${cpgtools_pyscript}" \
    -i ${cpgtools_dir}/input/${lab}_D6_vs_M8.bed \
    -g ${cpgtools_dir}/input/${lab}_D6_vs_M8.group \
    -o ${cpgtools_dir}/output/${lab}_D6_vs_M8.cpgtools &
  time python "${cpgtools_pyscript}" \
    -i ${cpgtools_dir}/input/${lab}_BC_vs_BL.bed \
    -g ${cpgtools_dir}/input/${lab}_BC_vs_BL.group \
    -o ${cpgtools_dir}/output/${lab}_BC_vs_BL.cpgtools & wait
done

labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
for lab in ${labs[@]}; do
  docker run --rm \
    -v /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/after/c_dmc/methylkit/input:/data/input \
    -v /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/after/c_dmc/methylkit/output:/data/output \
    -v /mnt/eqa/zhangyuanfeng/methylation/src/evaluation:/data/bin \
    methylkit:1.33.3 \
    Rscript /data/bin/methylkit.R -l ${lab}
done

labs=("BS1" "BS2" "BS3" "BS4" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1")
for lab in ${labs[@]}; do
  docker run --rm \
    -v /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/after/c_dmc/methylsig/input:/data/input \
    -v /mnt/eqa/zhangyuanfeng/methylation/best_pipeline/data/evaluated/after/c_dmc/methylsig/output:/data/output \
    -v /mnt/eqa/zhangyuanfeng/methylation/src/evaluation:/data/bin \
    methylsig:1.19.0 \
    Rscript /data/bin/methylsig.R -l ${lab}
done


# ! evaluate
p_types=("p" "q")
methylsig_methods=(
  "methylsig-beta-binomial"
  "methylsig-binomial"
  "methylsig-DSS"
)
treats=("before" "after")
for treat in "${treats[@]}"; do
for p_type in "${p_types[@]}"; do
  python src/evaluation/dmc_evaluate.py \
    -i best_pipeline/data/evaluated/${treat}/c_dmc/methylkit/output \
    -r quartet_reference/single_c/ensembl/dmc \
    -c best_pipeline/data/merged \
    -o best_pipeline/data/evaluated/${treat}/c_dmc/methylkit \
    -m methylkit -v info --p-type ${p_type} -g
  for methlysig_method in "${methylsig_methods[@]}"; do
    python src/evaluation/dmc_evaluate.py \
      -i best_pipeline/data/evaluated/${treat}/c_dmc/methylsig/output \
      -r quartet_reference/single_c/ensembl/dmc \
      -c best_pipeline/data/merged \
      -o best_pipeline/data/evaluated/${treat}/c_dmc/methylsig \
      -m ${methlysig_method} -v info --p-type ${p_type} -g
  done
  python src/evaluation/dmc_evaluate.py \
    -i best_pipeline/data/evaluated/${treat}/c_dmc/cpgtools/output \
    -r quartet_reference/single_c/ensembl/dmc \
    -c best_pipeline/data/merged \
    -o best_pipeline/data/evaluated/${treat}/c_dmc/cpgtools \
    -m cpgtools -v info --p-type ${p_type} -g
done
done