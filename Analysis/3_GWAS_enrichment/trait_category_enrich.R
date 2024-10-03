rm(list=ls())
library(tidyverse)
library(data.table)
library(parallel)

####################################################################
#
#             enrichment analysis by parent terms
#
###################################################################
########################### separate by tissue
cat <- fread('03_enrich_GWAS/gwas_catalog_trait-mappings_r2024-01-19.tsv')
names(cat) <- gsub(' ', '_', names(cat))
cat$Parent_term <- gsub(' ', '_', cat$Parent_term)
unique(cat$Parent_term)
cat <- filter(cat, Parent_term != 'NR')
table(cat$Parent_term)

GWAS.filter <- fread('./03_enrich_GWAS/GWAS_enrichment/GWAS_r2024-01-19.filter.txt')
GWAS <- fread('03_enrich_GWAS/gwas_catalog_v1.0.2-associations_e111_r2024-01-19.tsv', header = T, sep = "\t", na.strings = c("", "NA", "N/A"), quote="")


gwas.term <- lapply(unique(cat$Parent_term), function(term){
  
  EFO_id <- filter(cat, Parent_term == term)$EFO_URI
  GWAS.tmp <- filter(GWAS.filter, Maped_trait_URI %in% EFO_id)
  
  return(GWAS.tmp)
  
})
names(gwas.term) <- unique(cat$Parent_term)
saveRDS(gwas.term, './03_enrich_GWAS/GWAScatalog_201401/gwas.parent.term.rds')


## trait number in each parent term
gwas.term <- readRDS('./03_enrich_GWAS/GWAScatalog_201401/gwas.parent.term.rds')
num <- lapply(gwas.term, function(tmp){
  
  return(data.frame(
    trait.num = length(unique(tmp$Study_accession)),
    association.num = nrow(tmp)
  ))
  
  
})
num <- as.data.frame(do.call(rbind, num))
write.table(num, './03_enrich_GWAS/GWAScatalog_201401/GWAS_r2024-01-19.filter_parent_term_num.txt', quote = F, sep = '\t')



################## rename parent terms
gwas.term <- readRDS('./03_enrich_GWAS/GWAScatalog_201401/gwas.parent.term.rds')
gwas.term <-lapply(gwas.term, function(tmp){
  
  return(as.data.frame(tmp))
  
})

gwas.term <- as.data.frame(do.call(rbind, gwas.term)) %>% rownames_to_column(., var = 'Parent_term')
gwas.term$Parent_term <- gsub("\\.\\d+", "", gwas.term$Parent_term)

gwas.term$term <- gwas.term$Parent_term
gwas.term$term[gwas.term$term == 'Lipid_or_lipoprotein_measurement'] <- 'Lipid/lipoprotein markers'
gwas.term$term[gwas.term$term == 'Body_measurement'] <- 'Anthropometric traits'
gwas.term$term[gwas.term$term == 'Cardiovascular_measurement'] <- 'Cardiovascular indicators'
gwas.term$term[gwas.term$term == 'Hematological_measurement'] <- 'Hematological parameters'
gwas.term$term[grepl("Other", gwas.term$term)] <- 'Others'
gwas.term$term[gwas.term$term == 'Liver_enzyme_measurement'] <- 'Liver enzyme levels'
gwas.term$term[gwas.term$term == 'Biological_process'] <- "Physiological responses and habits"
gwas.term$term[gwas.term$term == 'Inflammatory_measurement'] <- "Inflammatory markers"
gwas.term$term <- gsub('_', ' ', gwas.term$term)


parent_term_idx <- unique(select(gwas.term, Parent_term, term))

a <- as.data.frame(gwas.term %>% 
                     group_by(term) %>% 
                     summarise(trait.num = n_distinct(Trait_mod), SNP.num = n_distinct(SNP)))

a <- merge(parent_term_idx, a, by = 'term')
write.table(a,'/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/03_enrich_GWAS/GWAScatalog_201401/parent_term_modify.txt', quote = F, sep = '\t', row.names = F)


### Count SNP number after parent term LD extend
outdir <- '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/03_enrich_GWAS/GWAScatalog_201401/'
term <- list.files(paste0(outdir, 'parent_term_mod_SNP_list/'))
num <- c()
for(x in term){
  
  ld <- fread(paste0(outdir, 'parent_term_mod_SNP_list_extend_LD0.5/', x, '.LD0.5.ld'))
  num <- data.frame(term = x, num = length(unique(ld$SNP_B[ld$R2 == 1]))) %>% rbind(num, .)
  
}

num$term <- gsub('_', ' ', num$term)
num$term[num$term == 'Lipid lipoprotein markers'] <- 'Lipid/lipoprotein markers'

num <- merge(a, num, by = 'term') 
unique <- unique(select(gwas.term, Parent_term, term))
unique
num <- merge(num, unique, by = 'term')
write.table(num, 'Data/parent_term_mod_SNP_list.txt', quote =F, sep = '\t', row.names = F)
###################################################################################



############################################### parent term LD extend ###############################################

outdir <- '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/03_enrich_GWAS/GWAScatalog_201401/'
lapply(unique(gwas.term$term), function(term.tmp){
  
  all <- filter(gwas.term, term == term.tmp) %>% select(., term, SNP)
  names(all) <- c('term', 'SNP')
  write.table(unique(all$SNP), paste0(outdir, 'parent_term_mod_SNP_list/', gsub(' ', '_', term.tmp) %>% gsub('/', '_', .)), row.names = F, sep = '\t', col.names =  F, quote = F)
  
})

## LD extend
dir <- '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/'
indir <- paste0(outdir, 'parent_term_mod_SNP_list/')
outdir1 <- paste0(outdir, 'parent_term_mod_SNP_list_extend_LD0.5/')

term <- list.files(paste0(outdir, 'parent_term_mod_SNP_list/'))

for(term.tmp in term){
  
  # run Plink to select SNP with LD > 0.5 within 500KB
  system(command =
           paste0('/workspace/rsrch1/ychen/miniconda3/bin/plink --bfile ', dir,
                  '1kg.v3/EUR --ld-snp-list ', indir, term.tmp, ' --r2 --ld-window 9999999 --ld-window-r2 0.5 --out ', outdir1,
                  term.tmp, '.LD0.5'))
}

################################################################################################################################### 

####################### Function enrichment 

rQTL.OR.enrichment <- function(term, id, qtl){
  
  all.ori <- fread(paste0('03_enrich_GWAS/GWAScatalog_201401/parent_term_mod_SNP_list_extend_LD0.5/', term, '.LD0.5.ld'))
  all <- filter(all.ori, R2 >= 1)
  gwas <- length(intersect(unique(id$rs_id_dbSNP151_GRCh38p7), unique(all$SNP_B)))
  bck <- length(unique(id$rs_id_dbSNP151_GRCh38p7))
  
  rqtl_gwas <- length(intersect(unique(qtl$rsid), unique(all$SNP_B)))
  rqtl_non_gwas <- length(unique(qtl$rsid)) - rqtl_gwas
  
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

### enrichment for different LD score
id <- fread('Data/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_maf.01.txt')
term.lsit <- list.files(paste0(outdir, 'parent_term_mod_SNP_list/'))

qtl <- fread(paste0('00_rQTL_mapping/04_Results_p/cis_rhyQTL_tissue/Liver.rQTL'))

all <- c()
for(term in term.list){

    tmp <- rQTL.OR.enrichment(term, id, qtl)
    all <- rbind(all, tmp)
}

write.table(all, './Result/03_enrich_GWAS/parent_term_mod_SNP_list_or/parent_term_enrich.txt', quote = F, row.names = F, sep = '\t')




############################################################# PLOT #############################################

theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8))


p1 <- ggplot(data = all, aes(x = odds.ratio, y = reorder(parent.term, odds.ratio))) + geom_point(size = 0.3, color= 'red') + theme +
  geom_vline(xintercept = 1) + labs(x = 'Enrichment', y = '') #+ xlim(0, 6)  
p1


library('grid')
pdf('Figure/03_enrich_GWAS/parent_term_liver_point.pdf', height = 2.6, width = 3.1)
print(p1, vp=viewport(.9, .9, x = .5, y = .5))
dev.off()

