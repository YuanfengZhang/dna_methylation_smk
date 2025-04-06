dbsnp_url="https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh38p7/VCF/00-common_all.vcf.gz"
dbsnp_target="resources/dbsnp/common_all.vcf.gz"
dbsnp_final="resources/dbsnp/common.dbsnp.vcf.gz"

grch38_gff_url="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/gencode.v47.basic.annotation.gff3.gz"
grch38_gff_target="resources/gencode.v47.basic.annotation.gff3.gz"

abismal_url="https://github.com/smithlabcode/abismal/releases/download/v3.2.4/abismal-3.2.4.tar.gz"
abismal_target="resources/abismal/bin/abismal"

mkdir -p resources/dbsnp

if [ ! -f "${dbsnp_final}" ]; then
    if ! aria2c -c -x 4 -k 1M "${dbsnp_url}" -o "${dbsnp_target}"; then
        echo "download failed: ${dbsnp_url}"
        exit 1
    fi
fi

zcat "${dbsnp_target}" |\
    awk 'BEGIN{OFS="\t"} $1 ~ /^#/ {print $0; next} {print "chr"$1,$2,$3,$4,$5,$6,$7,$8}' |\
    bgzip -@ 2 -i -I "${dbsnp_final}.tbi" >"${dbsnp_final}"

rm "${dbsnp_target}"


if [ ! -f "${grch38_gff_target}" ]; then
    if ! aria2c -c -x 4 -k 1M "${grch38_gff_url}" -o "${grch38_gff_target}"; then
        echo "download failed: ${grch38_gff_url}"
        exit 1
    fi
fi

mkdir -p resources/abismal

if [ ! -f "${abismal_target}" ]; then
    if ! aria2c -c -x 4 -k 1M "${abismal_url}" -o "resources/abismal-3.2.4.tar.gz"; then
        echo "download failed: ${abismal_url}"
        exit 1
    fi
    tar -xvzf "${abismal_target}.tar.gz"
    cd resources/abismal-3.2.4
    mkdir build && build
    ../configure --prefix ../baismal \
        CPPFLAGS="-I/opt/miniforge/envs/snakemake/include" \
        LDFLAGS="-L/opt/miniforge/envs/snakemake/lib"
    make -j 8 && make install
    rm "${abismal_target}.tar.gz"
    rm -rf resources/abismal-3.2.4
fi
