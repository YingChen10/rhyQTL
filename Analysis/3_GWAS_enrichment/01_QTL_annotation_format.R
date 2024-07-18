rm(list=ls())
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))

args <- commandArgs(TRUE)
tissue <- args[1]
chr <- args[2]

# chr <- 21
# tissue <- 'Liver'

geno_dir <- '/workspace/rsrch1/ychen/Projects/common_data/LDSC_resource/1000G_EUR_Phase3_plink/1000G.EUR.QC.'
bim <- fread(paste0(geno_dir, chr, ".bim"), head=F, stringsAsFactors=F, data.table=F)

rsid <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_maf.01.txt')

bim.1 <- bim[, c(1,4,2,3)]
bim.1$base <- 1
names(bim.1) <- c("CHR", "BP", "SNP", "CM", "base")
# annot <- cbind(bim.1, annot)

################################################################################
#
# observed QTLs
#
#################################################################################

# Import and format rQTL data
rqtl <- fread(paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/00_rQTL_mapping/04_Results_p/cis_rhyQTL_tissue/', tissue, '.rQTL')) %>% as.data.frame()

# top rQTLs
max <- data.frame(rqtl %>% group_by(gene) %>% filter(pval==min(pval)))
top <- max[which(!duplicated(max$gene)),]
top$tss_distance <- abs(top$pos - top$gene_start)
top$tss_distance[top$gene_sense == '-'] <- abs(top$gene_end[top$gene_sense == '-'] - top$pos[top$gene_sense == '-'])
rqtl_top <- select(top, tss_distance, MAF.Freq)

annot_rqtl <- top$rsid


# Import and format eQTL data
name_match <- fread("/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/03_enrich_GWAS/coloc/file_name.match.txt")
tissue_eqtl <- filter(name_match, rQTL_file == tissue)$eQTL_file
eqtl <- fread(paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/eQTLs/GTEx_Analysis_v8_eQTL/', tissue, '.v8.signif_variant_gene_pairs.txt.gz'))


# top eQTLs
max=data.frame(eqtl %>% group_by(gene_id) %>% filter(pval_nominal==min(pval_nominal)))
top=max[which(!duplicated(max$gene_id)),]
annot_eqtl_id <- top$variant_id
annot_eqtl <- filter(rsid, variant_id %in% annot_eqtl_id)$rs_id_dbSNP151_GRCh38p7


annot <- bim.1
annot$rQTL.obs <- 0
annot$rQTL.obs[match(annot_rqtl, annot$SNP)] <- 1

annot$eQTL.obs <- 0
annot$eQTL.obs[match(annot_eqtl, annot$SNP)] <- 1



indir <- '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/02_Functional_annotation_enrichment/00_input/ID_RDS/'
ann <- lapply(1:30, function(i){
  
  
  eQTL_random <- readRDS(paste0(indir, 'all_eQTL_random_ID/time.', i, '.rds'))
  rQTL_random <- readRDS(paste0(indir, 'all_rQTL_random_ID/time.', i, '.rds'))
  
  annot_eqtl_random <- filter(rsid, variant_id %in% eQTL_random[[tissue]])$rs_id_dbSNP151_GRCh38p7
  annot_rqtl_random <- filter(rsid, variant_id %in% rQTL_random[[tissue]])$rs_id_dbSNP151_GRCh38p7
  
  bim$rQTL_random_rand <- 0
  bim$rQTL_random_rand[match(annot_rqtl_random, bim$V2)] <- 1
  
  bim$eQTL_random_rand <- 0
  bim$eQTL_random_rand[match(annot_eqtl_random, bim$V2)] <- 1
  
  return(bim[,7:8])
  
})

result <- c()
for (i in 1:30) {
  result <- c(result, paste("rQTL.random.", i, sep = ""))
  result <- c(result, paste("eQTL.random.", i, sep = ""))
}

annot <- as.data.frame(do.call(cbind, ann))
names(annot) <- result


annot_dir <- paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/LDSC/01_annotation_top_pval/', tissue, '/')
if (!file.exists(annot_dir)) {dir.create(annot_dir)}

write.table(annot, paste0(annot_dir, 'erQTL.', chr, ".annot"), row = F, quo = F, sep = "\t")
system(paste0('gzip ', paste0(annot_dir, 'erQTL.', chr, ".annot")))