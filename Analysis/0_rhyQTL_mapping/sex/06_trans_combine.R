

rm(list=ls())

filter_regression <- function(male.reg, hanova.male){
  
  male.reg <- as.data.frame(do.call(rbind, male.reg)) %>% rownames_to_column(., var = 'idx')
  # male.reg$gene <- str_split_fixed(male.reg$idx, 'chr', 2)[ ,1]
  # male.reg$ID <- paste0("chr", str_split_fixed(male.reg$idx, 'chr', 2)[ ,2])
  
  male.reg$max_amp <- apply(select(male.reg, contains("amp.c")), 1, function(x) max(x))
  male.reg$max_amp_idx <- apply(select(male.reg, contains("amp.c_")), 1, function(x) which.max(x))
  male.reg$pval <- apply(select(male.reg, contains("pval_"), max_amp_idx), 1, function(x) (x[x[3]]))
  male.reg$qval <- apply(select(male.reg, contains("qval_"), max_amp_idx), 1, function(x) (x[x[3]]))
  
  male.reg <- filter(male.reg, pval < 5e-4)
  male.reg$idx.2 <- paste0(male.reg$idx, male.reg$max_amp_idx)
  
  
  
  hanova.male <- as.data.frame(do.call(rbind, hanova.male))[,c(1:4,8)] %>% rownames_to_column(., var = 'gene')
  hanova.male$gene <- paste0(str_split_fixed(hanova.male$gene, '\\.', 3)[ ,1], '.', str_split_fixed(hanova.male$gene, '\\.', 3)[ ,2])
  hanova.male$idx <- paste0(hanova.male$gene, '.', hanova.male$ID)
  
  male.reg <- merge(male.reg, hanova.male, by = 'idx')
  male.reg$idx.2 <- paste0(male.reg$idx, ':', male.reg$geno.1, ':', male.reg$geno.2)
  return(male.reg)
  
}

setwd('F:/Projects/Project03_human_circadian/rQTL/cis_QTL/11_Gender')

library(data.table)
library(tidyverse)

file <- list.files('05_compare_among_4_group')

all <- lapply(file, function(file.tmp) {
  
  message(file.tmp)
  
  male.reg <- readRDS(paste0('01_Rhythm_regression/Adipose-Visceral_Omentum.male/', file.tmp, '.rds'))
  hanova.male <- readRDS(paste0('04_hanova/Adipose-Visceral_Omentum.male/', file.tmp, '.rds'))
  # filter regerssion with p value < 5e-4
  male.reg.filter <- filter_regression(male.reg, hanova.male)
  
  
  female.reg <- readRDS(paste0('01_Rhythm_regression/Adipose-Visceral_Omentum.female/', file.tmp, '.rds'))
  hanova.female <- readRDS(paste0('04_hanova/Adipose-Visceral_Omentum.female/', file.tmp, '.rds'))
  # filter regerssion with p value < 5e-4
  female.reg.filter <- filter_regression(female.reg, hanova.male)
  
  # filter idx that have difference bwt 2 genotype group in both gender 
  intersect.idx.2 <- intersect(female.reg.filter$idx.2, male.reg.filter$idx.2)
  intersect.idx <- female.reg.filter$idx[female.reg.filter$idx.2 %in% intersect.idx.2]
  # commone in regression step
  length(intersect.idx)
  
  # dryR
  female.dryr <- fread(paste0('03_Rhythm_compare_multiple_times/Adipose-Visceral_Omentum.female/', file.tmp, '.filter')) %>% as.data.frame() %>% mutate(idx = paste0(gene, '.', ID))
  male.dryr <- fread(paste0('03_Rhythm_compare_multiple_times/Adipose-Visceral_Omentum.male.all/', file.tmp, '.filter')) %>% as.data.frame() %>% mutate(idx = paste0(gene, '.', ID))
  
  female.dryr <- filter(female.dryr, idx %in% intersect.idx) %>% mutate(gender = 'female')
  male.dryr <- filter(male.dryr, idx %in% intersect.idx) %>% mutate(gender = 'male')
  
  # common rhyQTL
  common.dryr <- intersect(female.dryr$idx, male.dryr$idx)
  common <- length(common.dryr)
  
  if(common == 0){return(NULL)}
  
  # gender specific idx 
  female.specific <- setdiff(unique(female.dryr$idx), common.dryr)
  male.specific <- setdiff(unique(male.dryr$idx), common.dryr)
  
  
  # test gender specific using 4 group comparing results
  compare.4 <- fread(paste0('05_compare_among_4_group/', file.tmp))
  compare.4$idx <- paste0(compare.4$gene, '.', compare.4$ID)
  female.specific.com.4 <- filter(compare.4, idx %in% female.specific) %>% filter(chosen_model == 28 | chosen_model == 33)
  male.specific.com.4 <- filter(compare.4, idx %in% male.specific) %>% filter(chosen_model == 18 | chosen_model == 23)
  
  count <- data.frame(file = file.tmp,
                      common = length(unique(filter(female.dryr, idx %in% common.dryr)$ID)),
                      female.specific = length(unique(female.specific.com.4$ID)),
                      male.specific = length(unique(male.specific.com.4$ID)),
                      female.specific.all = length(unique(filter(compare.4, idx %in% female.specific) %>% filter(chosen_model == 28 | chosen_model == 33) %>% .$ID)),
                      male.specific.all = length(unique(filter(compare.4, idx %in% male.specific) %>% filter(chosen_model == 18 | chosen_model == 23) %>% .$ID))) 
  
  
  rhyqtl <- rbind(data.frame(idx = common.dryr, chosen_model = 'same', chosen_model_BICW = 1, type = 'common'),
                  select(female.specific.com.4, idx, chosen_model, chosen_model_BICW) %>% mutate(type = 'female.specific'),
                  select(male.specific.com.4, idx, chosen_model, chosen_model_BICW) %>% mutate(type = 'male.specific'))
  
  
  
  result <- list(count = count,
                 rhyQTL = rhyqtl)
  
  return(result)
})

saveRDS(all, 'all.rds')

# count <- lapply(all, function(tmp){
#   
#   return(tmp$count)
#   
# })
# count <- as.data.frame(do.call(rbind, count))
# 
# all <- sum(count$common) + sum(count$female.specific) + sum(count$male.specific)
# sum(count$female.specific)/all
# sum(count$male.specific)/all


###################################### gender specific rhyQTL
library(tidyverse)
all <- readRDS('all.rds')
result <- lapply(all, function(tmp){
  
  return(tmp$rhyQTL)
  
})
result <- as.data.frame(do.call(rbind, result))
result$gene <- str_split_fixed(result$idx, '.chr', 2)[ ,1]
result$ID <- paste0('chr', str_split_fixed(result$idx, '.chr', 2)[ ,2])
head(result)

rqtl.count <- fread('../13_Results/00_rhyGene.rhyQTL.num.txt')

result.filter <- filter(result, chosen_model_BICW > 0.9)
annot <- fread('F:/Projects/Project03_human_circadian/rQTL/gencode.v26.GRCh38.genes.annot')
annot <- unique(annot[,6:7])
names(annot) <- c('gene','name')

result.filter <- merge(annot, result.filter, by = 'gene')
head(result.filter)
write.table(result.filter, 'gender_dif_rhyQTL.txt', quote = F, sep = '\t', row.names = F)

result.filter <- fread('gender_dif_rhyQTL.txt')
length(unique(result.filter$ID))

percentage <- table(result.filter$type)
pie(percentage)

count <- as.data.frame(table(result.filter$chosen_model))

sum(count$Freq)
# female specific
sum(count$Freq[count$Var1 == 28 | count$Var1 == 33])/sum(count$Freq)
# male specific
sum(count$Freq[count$Var1 == 18 | count$Var1 == 23])/sum(count$Freq)



