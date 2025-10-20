pairs=("D6_vs_D5" "D6_vs_F7" "D6_vs_M8" "BC_vs_BL")
for pair in "${pairs[@]}"; do
  python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_ref_prepare.py \
    pseudo_lab\
    -i /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/calibrated \
    -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE
done

labs=("BS0" "BS1" "BS2" "BS3" "BS4" "EM0" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1" "RM1")
pairs=("D6,D5" "D6,F7" "D6,M8" "BC,BL")
# tools=("cpgtools" "methylkit" "methylsig")
tool="methylsig"
for lab in "${labs[@]}"; do
  for pair in "${pairs[@]}"; do
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_prepare.py \
      -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
      -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/${tool}/input \
      -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/calibrated \
      -l ${lab} \
      -s ${pair} \
      -m ${tool} \
      -d 5 \
      -c chr17 \
      -v debug \
      -a prepare
  done
done

labs=("BS0" "BS1" "BS2" "BS3" "BS4" "EM0" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1" "RM1")
cpgtools_pyscript="/home/zhangyuanfeng/mambaforge/envs/cpgtools/bin/dmc_fisher.py"
cpgtools_dir="/mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/cpgtools"
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


labs=("BS0" "BS1" "BS2" "BS3" "BS4" "EM0" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1" "RM1")
for lab in ${labs[@]}; do
  docker run --rm \
    -v /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit/input:/data/input \
    -v /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit/output:/data/output \
    -v /mnt/eqa/zhangyuanfeng/methylation/src/evaluation:/data/bin \
    methylkit:1.33.3 \
    Rscript /data/bin/methylkit.R -l ${lab}
done

labs=("BS0" "BS1" "BS2" "BS3" "BS4" "EM0" "EM1" "EM2" "EM3" "EM4" "PS1" "PS2" "PS3" "RR1" "RM1")
for lab in ${labs[@]}; do
  docker run --rm \
    -v /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig/input:/data/input \
    -v /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig/output:/data/output \
    -v /mnt/eqa/zhangyuanfeng/methylation/src/evaluation:/data/bin \
    methylsig:1.19.0 \
    Rscript /data/bin/methylsig.R -l ${lab}
done


pairs=("D6_vs_D5" "D6_vs_F7" "D6_vs_M8" "BC_vs_BL")
p_types=("p" "q")
dmc_root_dir="/mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc"
cytosine_truset_dir="/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/calibrated"
dmc_truset_dir="/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc"

for pair in "${pairs[@]}"; do
  for p_type in "${p_types[@]}"; do
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_ref_prepare.py \
      merge_1040 \
      -m cpgtools --p-type ${p_type} \
      --bs0 "${dmc_root_dir}/cpgtools/output/BS0_${pair}.cpgtools.pval.txt" \
      --em0 "${dmc_root_dir}/cpgtools/output/EM0_${pair}.cpgtools.pval.txt" \
      -r "${cytosine_truset_dir}" \
      -o "${dmc_truset_dir}"
  done
done

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


pairs=("D6_vs_D5" "D6_vs_F7" "D6_vs_M8" "BC_vs_BL")
p_types=("p" "q")
dmc_root_dir="/mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc"
cytosine_truset_dir="/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/calibrated"
dmc_truset_dir="/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc"
for pair in "${pairs[@]}"; do
  for p_type in "${p_types[@]}"; do
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_ref_prepare.py \
      merge_1040 -m methylsig-beta-binomial --p-type ${p_type} \
      --bs0 "${dmc_root_dir}/methylsig/output/BS0_${pair}.beta_binomial.methylsig" \
      --em0 "${dmc_root_dir}/methylsig/output/EM0_${pair}.beta_binomial.methylsig" \
      -r "${cytosine_truset_dir}" \
      -o "${dmc_truset_dir}"
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_ref_prepare.py \
      merge_1040 -m methylsig-binomial --p-type ${p_type} \
      --bs0 "${dmc_root_dir}/methylsig/output/BS0_${pair}.binomial.methylsig" \
      --em0 "${dmc_root_dir}/methylsig/output/EM0_${pair}.binomial.methylsig" \
      -r "${cytosine_truset_dir}" \
      -o "${dmc_truset_dir}"
    python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_ref_prepare.py \
      merge_1040 -m methylsig-DSS --p-type ${p_type} \
      --bs0 "${dmc_root_dir}/methylsig/output/BS0_${pair}.dss.methylsig" \
      --em0 "${dmc_root_dir}/methylsig/output/EM0_${pair}.dss.methylsig" \
      -r "${cytosine_truset_dir}" \
      -o "${dmc_truset_dir}"
  done
done

dmc_root_dir="/mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc"
cytosine_truset_dir="/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/calibrated"
dmc_truset_dir="/mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc"

python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/cpgtools/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/cpgtools \
  -m cpgtools --p-type p -v debug
python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/cpgtools/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/cpgtools \
  -m cpgtools --p-type q -v debug
python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit \
  -m methylkit --p-type p -v debug
python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylkit \
  -m methylkit --p-type q -v debug
python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig \
  -m methylsig-beta-binomial --p-type p -v debug
python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig \
  -m methylsig-beta-binomial --p-type q -v debug
python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig \
  -m methylsig-binomial --p-type p -v debug
python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig \
  -m methylsig-binomial --p-type q -v debug
python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig \
  -m methylsig-DSS --p-type p -v debug
python /mnt/eqa/zhangyuanfeng/methylation/src/evaluation/dmc_evaluate.py \
  -i /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig/output \
  -r /mnt/eqa/zhangyuanfeng/methylation/quartet_reference/single_c/ensembl/dmc \
  -c /mnt/eqa/zhangyuanfeng/methylation/retrieve/CYTOSINE \
  -o /mnt/eqa/zhangyuanfeng/methylation/retrieve/dmc/methylsig \
  -m methylsig-DSS --p-type q -v debug

