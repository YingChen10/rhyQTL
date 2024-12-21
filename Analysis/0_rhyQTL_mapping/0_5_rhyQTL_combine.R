rm(list=ls())
library('tidyverse')
library('data.table')
library('parallel')

###################
#
# combine all split files
#
####################
setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/00_rQTL_mapping/')
tissue.list <- list.dirs('00_Genotype')[-1] %>% gsub('00_Genotype/', '', .)
file.list <- list.files('./00_Genotype/Adipose-Subcutaneous/') %>% gsub('.rds', '', .)


all <- mclapply(tissue.list, function(tissue){
  
  message(tissue)
  file.tmp <- lapply(file.list, function(file){
    
    # regerssion
    regression.tmp <- readRDS(paste0('01_Rhythm_regression/', tissue, '/', file, '.rds'))
    
    regression.tmp <- lapply(regression.tmp, function(tmp){
      
      tmp <- rownames_to_column(tmp, var = 'ID')
      return(tmp)
      
    })
    
    regression.tmp <- as.data.frame(do.call(rbind,  regression.tmp)) %>% rownames_to_column(., var = 'gene')
    regression.tmp$gene <- paste0(str_split_fixed(regression.tmp$gene, '\\.', 3)[ ,1], '.', str_split_fixed(regression.tmp$gene, '\\.', 3)[ ,2])
    regression.tmp$max_amp <- apply(select(regression.tmp, contains("amp.c")), 1, function(x) max(x))
    regression.tmp$max_amp_idx <- apply(select(regression.tmp, contains("amp.c_")), 1, function(x) which.max(x))
    regression.tmp$pval <- apply(select(regression.tmp, contains("pval_"), max_amp_idx), 1, function(x) (x[x[3]]))
    regression.tmp$qval <- apply(select(regression.tmp, contains("qval_"), max_amp_idx), 1, function(x) (x[x[3]]))
    regression.tmp <- filter(regression.tmp, pval < 5e-4)
    
    # hanaove
    hanova.tmp <- readRDS(paste0('04_hanova/', tissue, '/', file, '.rds'))
    hanova.tmp <- as.data.frame(do.call(rbind,  hanova.tmp))
    hanova.tmp <- hanova.tmp[,c(1:4, 8)]
    hanova.tmp <- rownames_to_column(hanova.tmp, var = 'gene')
    hanova.tmp$gene <- paste0(str_split_fixed(hanova.tmp$gene, '\\.', 3)[ ,1], '.', str_split_fixed(hanova.tmp$gene, '\\.', 3)[ ,2])
    
    # combine
    combine <- merge(regression.tmp, hanova.tmp, by = c('gene', 'ID'))
    
    
    #compare 
    compare <- fread(paste0('03_Rhythm_compare_multiple_times/', tissue, '/', file))
    combine <- merge(combine, compare,  by = c('gene', 'ID'))
    
    return(combine)
    
  })
  
  file.tmp <- as.data.frame(do.call(rbind, file.tmp)) 
  file.tmp <- filter(file.tmp, gtest.p.value < 0.05)
  file.tmp$p_adj <- p.adjust(file.tmp$HANOVA.norm, method = "BH")
  file.tmp <- filter(file.tmp, p_adj < 0.01)
  
  file.tmp$tissue <- tissue
  
  return(file.tmp)
  
}, mc.cores = 5)
names(all) <- tissue.list
saveRDS(all, '../13_Results/rhyQTL_rhyGene.rds') 

################################################
#
# combine all rhyQTLs and rhyGenes
#
################################################

dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/"
cat <- "cis_QTL/"


annot <- read.table(paste0(dir, 'gencode.v26.GRCh38.genes.annot'))
annot <- unique(annot[, c('V6', "V7")])
names(annot) <- c('gene', 'name')

outdir <- paste0(dir, cat, "13_Results/")
if (!file.exists(outdir)) {dir.create(outdir, recursive = TRUE)}

outdir.all <- paste0(dir, cat, "13_Results/00_all_rhyQTL/")
if (!file.exists(outdir.all)) {dir.create(outdir.all, recursive = TRUE)}

outdir.lead <- paste0(dir, cat, "13_Results/00_lead_rhyQTL/")
if (!file.exists(outdir.lead)) {dir.create(outdir.lead, recursive = TRUE)}

all <- readRDS('../13_Results/rhyQTL_rhyGene.rds') 
combine <- data.frame()
for(i in 1:length(all)){
  
  tissue.list <- all[[i]]
  tissue.list <- merge(annot, tissue.list, by = 'gene')
  
  write.table(tissue.list, paste0(outdir.all, names(all)[i], '.txt'), quote = F, sep = '\t', row.names = F)
  
  
  result <- tissue.list %>%
    group_by(gene) %>%
    slice_min(order_by = HANOVA.norm, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  write.table(result, paste0(outdir.lead, names(all)[i], '.txt'), quote = F, sep = '\t', row.names = F)
  
  combine <- unique(select(result, gene, name)) %>% mutate(tissue = names(all)[i]) %>% rbind(combine, .)
  
  
}

write.table(combine, paste0(outdir,'00_tissue_rhyGene.name.list'), quote = F, sep = '\t', row.names = F)
