
# generate motif file as the input of meme 
/opt/meme-5.4.1/scripts/chen2meme homer_motif.txt > homer_motif_meme.txt

for chr in {1..22}
do


samtools faidx GRCh38.primary_assembly.genome.fa chr${chr} > ${chr}.fa
# scan on genome
fimo -o fimo_out_${chr} homer_motif_meme.txt ${chr}.fa
# format
awk '{print $2,$3,$4,$1,$5,$6,$7,$8,$9}' fimo_out_${chr}/fimo.txt|grep -v '#' |grep -v "motif_id"| sed "s/\s/\t/g"| sort -k1,1 -k2,2n > fimo_out_${chr}/fimo.bed
# intersect
bedtools intersect -a GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_maf.01.bed -b fimo_out_${chr}/fimo.bed -wa -wb > fimo_out_${chr}/snp.bed

done

Rscript TF_motif_enrichment.R
Rscript Figure_TF_motif_enrich.R
