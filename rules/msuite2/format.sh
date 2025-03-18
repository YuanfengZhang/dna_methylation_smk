cd $1

echo "rename bedgraph"
mv Msuite2.CpG.meth.bedgraph.gz "$2.msuite2.bedgraph.gz"

echo "reheader bam files"
bam_files=("Msuite2.final" "Lambda.srt" "pUC19.srt")
for bam_f in ${bam_files[@]}; do
    samtools view -H "${bam_f}.bam" |\
    sed -e 's/chrL/lambda/' -e 's/chrP/pUC19/' > "${bam_f}.header.sam"
    samtools reheader "${bam_f}.header.sam" "${bam_f}.bam" > "${bam_f}.reheader.bam"
    rm "${bam_f}.bam" "${bam_f}.bam.bai"
done

mv Lambda.srt.reheader.bam lambda.bam
mv pUC19.srt.reheader.bam pUC19.bam

echo "merge"
samtools merge -o "$2.bam" \
    -f -p -@ $3 \
    Msuite2.final.reheader.bam \
    lambda.bam \
    pUC19.bam

samtools index -@ $3 "$2.bam"
samtools index -@ $3 lambda.bam
samtools index -@ $3 pUC19.bam

rm Msuite2.final.reheader.bam
