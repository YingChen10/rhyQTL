dir="/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/"
ldsc="/workspace/rsrch1/ychen/Projects/common_data/LDSC_resource"

tissue=$1
trait=$2

mkdir -p ${dir}/Result/03_enrich_GWAS/LDSC/02_Heritability2/${tissue}

ldsc_env="/workspace/rsrch1/ychen/miniconda3/envs/ldsc/bin/python"
${ldsc_env} /workspace/rsrch1/ychen/bin/ldsc/ldsc.py \
	--h2 ${dir}/03_enrich_GWAS/LDSC/sumstats/GWAScata_sumstats/${trait}.sumstats.gz \
	--ref-ld-chr ${dir}/Result/03_enrich_GWAS/LDSC/01_annotation/${tissue}/erQTL. \
	--w-ld-chr ${ldsc}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
	--overlap-annot \
	--frqfile-chr ${ldsc}/1000G_Phase3_frq/1000G.EUR.QC. \
	--out ${dir}/Result/03_enrich_GWAS/LDSC/02_Heritability2/${tissue}/${trait}.ldsc \
	--print-coefficients --print-delete-vals
