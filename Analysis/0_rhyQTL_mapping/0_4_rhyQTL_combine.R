rm(list = ls()) 
library('tidyverse')
library(parallel)
library('data.table')


################################ 1. 03_Rhythm_compare_multiple_times results combine ############################################
setwd <- '/workspace/rsrch1/ychen/Projects/rQTL/cis_rQTL'
file <- list.files('/workspace/rsrch1/ychen/Projects/rQTL/GTEx_nor_expression/', pattern = 'txt')
tissue <- gsub('.txt', '', file)
tissue <- tissue[tissue != 'Kidney-Cortex']

rQTL <- mclapply(tissue, function(x){
  
  all <- c()
  file <- list.files(paste0('./00_rQTL_mapping/03_Rhythm_compare_multiple_times/', x), pattern = 'filter')
  for(y in file){
    
    # rhymthmic compare results
    tmp <- read.table(paste0('./00_rQTL_mapping/03_Rhythm_compare_multiple_times/', x, '/', y), header = T)
    all <- rbind(all, tmp)
  }
  
  return(all)
  
}, mc.cores = 45)

names(rQTL) <- tissue
saveRDS(rQTL, file = "00_rQTL_mapping/04_Result/03_compare_multiple_times_rQTL.rds")
rQTL <- readRDS("00_rQTL_mapping/04_Result/03_compare_multiple_times_rQTL.rds")
###########################################################################################################



############################################ 2. 01_Rhythm_regression results combine #################################################
split <- read.table('./split_pos.file')

regression <- mclapply(tissue, function(x){
  
  all <- c()
 
  for(y in split$V1){
    
    # rhymthmic compare results
    tmp <- readRDS(paste0('./00_rQTL_mapping/01_Rhythm_regression/', x, '/', y, '.rds'))
    
    tmp.1 <- lapply(tmp, function(gene.list){
      
      gene.list <- rownames_to_column(gene.list, var = 'ID')
      
    })
    
    tmp.1 <- as.data.frame(do.call(rbind, tmp.1))
    tmp.1 <- rownames_to_column(tmp.1, var = 'gene')
    tmp.1$gene <- paste0(str_split_fixed(tmp.1$gene, '\\.', 3)[ ,1], '.', str_split_fixed(tmp.1$gene, '\\.', 3)[ ,2])
    
    all <- rbind(all, tmp.1)
  }
  
  return(all)
  
}, mc.cores = 45)

names(regression) <- tissue
saveRDS(regression, file = "00_rQTL_mapping/04_Result/01_rhythm_regression.rds")
####################################################

############################################# 3. merge step 1 and 2 #######################
results <- mclapply(tissue, function(x){
  
  all <- merge(regression[[x]], rQTL[[x]], by = c('gene', 'ID'))
  return(all)
  
}, mc.cores = 45)
names(results) <- tissue
saveRDS(results, file = "00_rQTL_mapping/04_Result/03_compare_multiple_times_rQTL_pvalue.rds")


## retrieve all cis-rQTL, use p value cutoff 10e-4
all.pval <- mclapply(results, function(tissue.list){
  
  tissue.list$max_amp <- apply(select(tissue.list, contains("amp.c")), 1, function(x) max(x))
  tissue.list$max_amp_idx <- apply(select(tissue.list, contains("amp.c_")), 1, function(x) which.max(x))
  tissue.list$pval <- apply(select(tissue.list, contains("pval_"), max_amp_idx), 1, function(x) (x[x[3]]))
  tissue.list$qval <- apply(select(tissue.list, contains("qval_"), max_amp_idx), 1, function(x) (x[x[3]]))
  tissue.list <- filter(tissue.list, pval <= 0.0001)
  
  return(tissue.list)
  
})

saveRDS(all.pval, file = "00_rQTL_mapping/04_Results/cis_rQTL.rds")
