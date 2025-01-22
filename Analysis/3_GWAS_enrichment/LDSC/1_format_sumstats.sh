#conda activate
#source activate ldsc
dir="/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/03_enrich_GWAS/LDSC/sumstats/"

python /workspace/rsrch1/ychen/bin/ldsc/munge_sumstats.py \
--sumstats /mnt/disk5_7T/LDSC/sumstats/Lipid_or_lipoprotein_measurement.GCST002221.Cholesterol-total.sumstats \
--merge-alleles /workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/03_enrich_GWAS/LDSC/w_hm3.snplist \
--out $dir/Lipid_or_lipoprotein_measurement.GCST00221.Cholesterol-total \
--a1-inc --chunksize 500000
