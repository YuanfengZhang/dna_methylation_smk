dbsnp_url="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz"

if command -v aria2c >/dev/null 2>&1; then
  aria2c \
      -c -x 4 -k 1M \
      "${dbsnp_url}" \
      -o "resources/common_all.vcf.gz"
else
  wget -c "${dbsnp_url}" -o "resources/common_all.vcf.gz"
fi

cd resources
zcat common_all.vcf.gz |\
    awk 'BEGIN{OFS="\t"} $1 ~ /^#/ {print $0; next} {print "chr"$1,$2,$3,$4,$5,$6,$7,$8}' |\
    bgzip -@ 2 -i -I common.dbsnp.vcf.gz.tbi >common.dbsnp.vcf.gz
