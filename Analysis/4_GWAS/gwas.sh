
dir="/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/06_VIVA/GWAS/"

### Step 0 ### 
cd 0_Input_format
# format VCF and phenotype file
Rscript input_format.R
Rscript combine_array_snp.R


plink --vcf array.snv.pheno.geno.matched.vcf --make-bed --out array.snv.pheno.geno.matched
## generate .fam file 
Rscript input_format.R

# After format the genotype and phenotype data we perfom GWAS following https://github.com/MareesAT/GWA_tutorial
