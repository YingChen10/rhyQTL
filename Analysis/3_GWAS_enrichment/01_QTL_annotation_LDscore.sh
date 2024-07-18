dir="/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/"
ldsc="/workspace/rsrch1/ychen/Projects/common_data/LDSC_resource/"
ldsc_env="/workspace/rsrch1/ychen/miniconda3/envs/ldsc/bin/python"

tissue='Liver'
for chr in {1..22}
do
  #############################################
  # 0. Format annotation for QTL
  #############################################
  # conda deactivate
  Rscript /workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Script/03_enrich_GWAS/2024-01/20240430_LDSC/01_QTL_annotation.R ${tissue} ${chr}

  #############################################
  # 1. Computing LD scores with an annot file
  #############################################
  ${ldsc_env} /workspace/rsrch1/ychen/bin/ldsc/ldsc.py \
        --l2 \
        --bfile ${ldsc}/1000G_EUR_Phase3_plink/1000G.EUR.QC.${chr} \
        --ld-wind-cm 1 \
        --annot ${dir}/Result/03_enrich_GWAS/LDSC/01_annotation/${tissue}/erQTL.${chr}.annot.gz \
        --out ${dir}/Result/03_enrich_GWAS/LDSC/01_annotation/${tissue}/erQTL.${chr} \
        --print-snps ${ldsc}/hapmap3_snps/hm.${chr}.snp  
done


