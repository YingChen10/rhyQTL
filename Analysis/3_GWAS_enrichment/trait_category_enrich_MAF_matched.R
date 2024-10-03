library(tidyverse)
library(data.table)

setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')

or.bin.calculate <- function(term, MAF.Freq, qtl, id){
  
  all.ori <- fread(paste0('03_enrich_GWAS/GWAScatalog_201401/parent_term_mod_SNP_list_extend_LD0.5/', term, '.LD0.5.ld'))
  all <- filter(all.ori, R2 >= 1)
  
  or <- c()
  for(i in 2:5){
    
  MAF.Freq.min <- MAF.Freq[i-1]
  MAF.Freq.max <- MAF.Freq[i]
    
  tmp <- filter(qtl, MAF.Freq > MAF.Freq.min & MAF.Freq <= MAF.Freq.max)
  bck.tmp <- filter(id, MAF.Freq > MAF.Freq.min & MAF.Freq <= MAF.Freq.max)

  rqtl_gwas <- length(intersect(unique(tmp$rsid), unique(all$SNP_B)))
  
  # ranodmly choose MAF mateched SNPs
  maf.mateched.snp <- sample(bck.tmp$rs_id_dbSNP151_GRCh38p7, length(unique(tmp$ID)))
  random_snp_gwas <- length(intersect(unique(maf.mateched.snp), unique(all$SNP_B)))
  
  or <- data.frame(parent.term = term,
                   rqtl.gwas =  rqtl_gwas, 
                   random.snp.gwas = random_snp_gwas,
                   MAF.bin = paste0(MAF.Freq.min, '-', MAF.Freq.max)) %>% rbind(or, .)
  
  }
  
  return(or)
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

sum <- c()
for(i in 1:30){
  
  print(i)
  results <- lapply(term.list, function(term){
    
    or.bin.calculate(term, MAF.Freq, qtl, id)
    
  })
  
  results <- as.data.frame(do.call(rbind, results))
  sum <- results %>% group_by(parent.term) %>%
    summarise(rqtl.gwas.num = sum(rqtl.gwas), random.snp.gwas.num = sum(random.snp.gwas)) %>% 
    mutate(time = i) %>%
    rbind(sum, .)
  
}

sum$ratio <- sum$rqtl.gwas.num/sum$random.snp.gwas.num
write.table(sum, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/03_enrich_GWAS/GWAScatalog_201401/ratio_random_snp.txt', quote = F, sep = '\t', row.names = F)

ratio <- sum %>% group_by(parent.term) %>% summarise(mean = mean(ratio), median = median(ratio), sd = sd(ratio))
ratio <- filter(ratio, parent.term != 'Others')
ratio$parent.term <- gsub('_', ' ', ratio$parent.term)
theme <- theme_bw() +
  theme(axis.line = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) +
  theme(legend.position = "none")



# ggplot(data = ratio, aes(x = mean, y = reorder(parent.term, mean))) +
#   geom_col(fill = 'red', width = 0.5) + 
#   geom_errorbar(aes(xmin = mean - sd, xmax = mean + sd), width = 0) + 
#   theme + 
#   labs(x = 'Enrichment', y = '') +
#   geom_vline(xintercept = 1, linetype = "dashed", color = 'grey20')

p1 <- ggplot(data = ratio, aes(x = median, y = reorder(parent.term, median))) +
  geom_col(fill = 'red', width = 0.5) + 
  geom_errorbar(aes(xmin = mean - sd, xmax = mean + sd), width = 0) + 
  theme + 
  labs(x = 'Enrichment', y = '') +
  geom_vline(xintercept = 1, linetype = "dashed", color = 'grey20')
         

pdf("Figures/03_GWAS_enrichment/GWAScatalog_201401/gwas_enrichment_MAF_match_random_snp.pdf", useDingbats = F, width = 14, height = 3)
print(p1, vp=viewport(.22, 1, x = 0.2, y = 0.5))
dev.off()
