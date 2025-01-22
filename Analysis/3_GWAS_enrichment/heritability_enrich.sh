#####################################
# 0.format GWAS sumstats
#####################################

trait_ids=(
    "GCST000755"
"GCST002216"
"GCST002221"
"GCST90092808"
"GCST90092809"
"GCST90092875"
"GCST90092928"
"GCST90092975"
"GCST90092980"
"GCST90092981"
"GCST90092987"
"GCST90093002"
"GCST90267273"
"GCST90269582"
"GCST90269553"
)


for trait_id in "${trait_ids[@]}"
do
    Rscript 00.0_format_sumstats.R "${trait_id}"
done

traits=(
    "Lipid_or_lipoprotein_measurement.GCST000755.HDL-cholesterol"
"Lipid_or_lipoprotein_measurement.GCST002216.Triglycerides"
"Lipid_or_lipoprotein_measurement.GCST002221.Cholesterol-total"
"Lipid_or_lipoprotein_measurement.GCST90092808.Apolipoprotein-A1-levels"
"Lipid_or_lipoprotein_measurement.GCST90092809.Apolipoprotein-B-levels"
"Lipid_or_lipoprotein_measurement.GCST90092875.Concentration-of-large-VLDL-particles"
"Lipid_or_lipoprotein_measurement.GCST90092928.Monounsaturated-fatty-acid-levels"
"Lipid_or_lipoprotein_measurement.GCST90092975.Concentration-of-small-VLDL-particles"
"Lipid_or_lipoprotein_measurement.GCST90092980.Saturated-fatty-acid-levels"
"Lipid_or_lipoprotein_measurement.GCST90092981.Ratio-of-saturated-fatty-acids-to-total-fatty-acids"
"Lipid_or_lipoprotein_measurement.GCST90092987.Total-fatty-acid-levels"
"Lipid_or_lipoprotein_measurement.GCST90093002.Average-diameter-for-VLDL-particles"
"Lipid_or_lipoprotein_measurement.GCST90267273.LDL-weighted-GWA"
"Lipid_or_lipoprotein_measurement.GCST90269553.Linoleic-acid-to-total-fatty-acids-percentage-UKB-data-field-23456"
"Lipid_or_lipoprotein_measurement.GCST90269582.Cholesteryl-esters-in-chylomicrons-and-extremely-large-VLDL-UKB-data-field-23485"
)

dir="/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/"

for trait in "${traits[@]}"
do

python /workspace/rsrch1/ychen/bin/ldsc/munge_sumstats.py --sumstats /mnt/disk/LDSC/sumstats/${trait}.sumstats --merge-alleles ${dir}/03_enrich_GWAS/LDSC/w_hm3.snplist --out ${dir}/03_enrich_GWAS/LDSC/sumstats/GWAScata_sumstats/${trait} --a1-inc --chunksize 500000

done

#####################################
# 1. calculate LD scores
#####################################

ldsc="/workspace/rsrch1/ychen/Projects/common_data/LDSC_resource/"
ldsc_env="/workspace/rsrch1/ychen/miniconda3/envs/ldsc/bin/python"

tissue='Liver'
for chr in {1..22}
do
  
  # 0. Format annotation for QTL
  # conda deactivate
  Rscript ${dir}/Script/03_enrich_GWAS/2024-01/20240430_LDSC/01_QTL_annotation_format.R ${tissue} ${chr}

  # 1. Computing LD scores with an annot file
  ${ldsc_env} /workspace/rsrch1/ychen/bin/ldsc/ldsc.py \
        --l2 \
        --bfile ${ldsc}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
        --ld-wind-cm 1 \
        --annot ${dir}/Result/03_enrich_GWAS/LDSC/01_annotation/${tissue}/erQTL.${chr}.annot.gz \
        --out ${dir}/Result/03_enrich_GWAS/LDSC/01_annotation/${tissue}/erQTL.${chr} \
        --print-snps ${ldsc}/hapmap3_snps/hm.${chr}.snp  
done


#############################################
# 2. Heritability
#############################################
tissue='Liver'

for trait in "${traits[@]}"
do
    sh 02_Heritability.sh "${tissue}" "${trait}"
done



