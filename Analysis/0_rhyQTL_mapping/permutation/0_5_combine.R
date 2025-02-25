rm(list=ls())
library('tidyverse')
library('data.table')
library('parallel')

###################
#
# combine all split files
#
####################
setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/12_permutation/')
tissue.list <- list.dirs('00_Genotype')[-1] %>% gsub('00_Genotype/', '', .)
file.list <- list.files('./00_Genotype/Adipose-Subcutaneous/') %>% gsub('.rds', '', .)


dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/"
cat <- "cis_QTL/12_permutation/"
outdir <- paste0(dir, cat, "05_combine/")
if (!file.exists(outdir)) {dir.create(outdir, recursive = TRUE)}

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
    #regression.tmp <- filter(regression.tmp, pval < 5e-4)

    # hanova
    hanova.tmp <- readRDS(paste0('04_hanova/', tissue, '/', file, '.rds'))
    hanova.tmp <- as.data.frame(do.call(rbind,  hanova.tmp))
    hanova.tmp <- hanova.tmp[,c(1:4, 8)]
    hanova.tmp <- rownames_to_column(hanova.tmp, var = 'gene')
    hanova.tmp$gene <- paste0(str_split_fixed(hanova.tmp$gene, '\\.', 3)[ ,1], '.', str_split_fixed(hanova.tmp$gene, '\\.', 3)[ ,2])

    # combine
    combine <- merge(regression.tmp, hanova.tmp, by = c('gene', 'ID'))


    # compare
    compare <- fread(paste0('03_Rhythm_compare_multiple_times/', tissue, '/', file))
    combine <- merge(combine, compare,  by = c('gene', 'ID'))

    return(combine)

  })

  file.tmp <- as.data.frame(do.call(rbind, file.tmp))

  file.tmp <- filter(file.tmp, pval <= 5e-4)
  file.tmp <- filter(file.tmp, gtest.p.value < 0.05)

  # length(unique(file.tmp$gene))
  # length(unique(file.tmp$ID))

  file.tmp$p_adj <- p.adjust(file.tmp$HANOVA.norm, method = "BH")
  file.tmp <- filter(file.tmp, p_adj <= 0.05)

  # length(unique(file.tmp$gene))
  # length(unique(file.tmp$ID))

  file.tmp$tissue <- tissue

  return(file.tmp)

}, mc.cores = 25)
names(all) <- tissue.list
saveRDS(all, paste0(outdir, 'rhyQTL_rhyGene.rds'))

################################################
#
# combine all rhyQTLs and rhyGenes
#
################################################


all <- readRDS(paste0(outdir, 'rhyQTL_rhyGene.rds'))
annot <- read.table('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/all_gene.annot')
annot <- unique(annot[, c('V6', "V7")])
names(annot) <- c('gene', 'name')

outdir.all <- paste0(outdir, "00_all_rhyQTL/")
if (!file.exists(outdir.all)) {dir.create(outdir.all, recursive = TRUE)}

outdir.lead <- paste0(outdir, "00_lead_rhyQTL/")
if (!file.exists(outdir.lead)) {dir.create(outdir.lead, recursive = TRUE)}

for(i in 1:length(all)){

  tissue.list <- all[[i]]
  tissue.list <- merge(annot, tissue.list, by = 'gene')

  write.table(tissue.list, paste0(outdir.all, names(all)[i], '.txt'), quote = F, sep = '\t', row.names = F)

}
