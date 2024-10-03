
library(tidyverse)
library(data.table)

setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')

or.bin.calculate <- function(qtl, MAF.Freq.min, MAF.Freq.max, id, term){
  
  tmp <- filter(qtl, MAF.Freq > MAF.Freq.min & MAF.Freq <= MAF.Freq.max)
  bck.tmp <- filter(id, MAF.Freq > MAF.Freq.min & MAF.Freq <= MAF.Freq.max)
  
  all.ori <- fread(paste0('03_enrich_GWAS/GWAScatalog_201401/parent_term_mod_SNP_list_extend_LD0.5/', term, '.LD0.5.ld'))
  all <- filter(all.ori, R2 >= 1)
  
  gwas <- length(intersect(unique(bck.tmp$rs_id_dbSNP151_GRCh38p7), unique(all$SNP_B)))
  bck <- length(unique(bck.tmp$rs_id_dbSNP151_GRCh38p7))
  
  rqtl_gwas <- length(intersect(unique(tmp$rsid), unique(all$SNP_B)))
  rqtl_non_gwas <- length(unique(tmp$rsid)) - rqtl_gwas
  
  non_rqtl_gwas <- gwas - rqtl_gwas
  non_rqtl_non_gwas <- bck - gwas - rqtl_non_gwas
  
  num <- c(rqtl_gwas, rqtl_non_gwas, non_rqtl_gwas, non_rqtl_non_gwas)
  
  data <- matrix(num, nrow = 2) 
  # colnames(data) <- c("rqtl", "non_rqtl") 
  # rownames(data) <- c("motif", "non_motif")
  
  OR <- data.frame(parent.term = term,
                   rqtl.gwas =  rqtl_gwas, 
                   rqtl.non.gwas = rqtl_non_gwas,  
                   non.rqtl.gwas = non_rqtl_gwas, 
                   non.rqtl.non.gwas = non_rqtl_non_gwas,
                   odds.ratio = as.numeric(fisher.test(data)$estimate),
                   p.value = as.numeric(fisher.test(data)$p.value))
  
  return(OR)
}

################## all bck snp
id.1 <- fread('Data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_maf.01.txt')
id <- fread('Data/02_GTEx_MAF.59_SNP_MAF_TSS_distance.txt')
id <- merge(id, id.1, by.x = 'ID', by.y = 'variant_id')


length(intersect(id$rs_id_dbSNP151_GRCh38p7, qtl$rsid))
length(unique(qtl$rsid))

################### liver rhyQTL snp
qtl <- fread(paste0('00_rQTL_mapping/04_Results_p/cis_rhyQTL_tissue/Liver.rQTL'))
range(qtl$MAF.Freq)
length(qtl$ID[qtl$MAF.Freq < 0.1]) # only 6 SNPs with MAF < 0.1, they will be assigned to the bin: 0.1 < MAF <= 0.2
qtl$MAF.Freq[qtl$MAF.Freq < 0.1] <- 0.1

outdir <- '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/03_enrich_GWAS/GWAScatalog_201401/'
term.list <- list.files(paste0(outdir, 'parent_term_mod_SNP_list/'))

#################### set MAF bin
MAF.Freq <- c(0.1, 0.2, 0.3, 0.4, 0.5)
results <- c()
for(i in 2:5){
  
  print(i)
  MAF.Freq.min <- MAF.Freq[i-1]
  MAF.Freq.max <- MAF.Freq[i]
  
  all <- c()
  for(x in term.list){
    
    or <- or.bin.calculate(qtl, MAF.Freq.min, MAF.Freq.max, id, x)
    all <- rbind(all, or)
    
  }
  all$MAF.Freq.bin <- paste0(MAF.Freq.min, '-', MAF.Freq.max)
  
  results <- rbind(results, all)
  
}

results <- filter(results, parent.term != 'Others')
write.table(results, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/03_enrich_GWAS/GWAScatalog_201401/MAF_bin_odds_ratio.txt', quote = F, sep = '\t', row.names = F)

unique(results$MAF.Freq.bin)

results$parent.term <- gsub('_', ' ', results$parent.term)

theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) +
  theme(legend.position = "none")


for(i in 1:length(unique(results$MAF.Freq.bin))){
  
  assign(paste0('p', i), 
         ggplot(data = filter(results, MAF.Freq.bin == unique(results$MAF.Freq.bin)[i]), aes(x = odds.ratio, y = reorder(parent.term, odds.ratio))) +
           geom_point(size = 0.3, color = 'red') + theme + labs(x = 'Enrichment', y = '', title = unique(results$MAF.Freq.bin)[i]) +
    geom_vline(xintercept = 1, linetype = "dashed", color = 'grey20'))
       
}

pdf("Figures/03_GWAS_enrichment/GWAScatalog_201401/gwas_enrichment_MAF_bin.pdf", useDingbats = F, width = 14, height = 3)
print(p1, vp=viewport(.24, 1, x = 0.1, y = 0.5))
print(p2, vp=viewport(.24, 1, x = 0.35, y = 0.5))
print(p3, vp=viewport(.24, 1, x = 0.6, y = 0.5))
print(p4, vp=viewport(.24, 1, x = 0.85, y = 0.5))
dev.off()
