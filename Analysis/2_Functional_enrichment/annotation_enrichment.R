## Function get control SNP count
Get.Count <- function(baseline, feature){
  
  ### baseline
  # add tissue information
  baseline <- baseline %>% mutate(tissue = gsub('\\.[0-9]+', '', row.names(.)))
  # calculate total number of QTLs in baseline
  baseline <- as.data.frame(baseline %>% group_by(tissue) %>% summarise(sum = sum(num)))
  
  
  ### feature
  feature <- feature %>% mutate(tissue = gsub('\\.[0-9]+', '', row.names(.)))
  
  count <- gather(feature, Feature, Count, enhancer_d:synonymous_variant_d) 
  
  # merge baseline and feature
  merge <- merge(count, baseline, by = 'tissue')
  names(merge) <- c('tissue', 'Feature', 'Count', 'Baseline')
  
  return(merge)
}

# Function get observed SNP count
Get.obs.Count <- function(baseline, feature){
  
  ### baseline
  # add tissue information
  baseline <- baseline %>% mutate(tissue = gsub('\\.[0-9]+', '', row.names(.)))
  # calculate total number of QTLs in baseline
  baseline <- as.data.frame(filter(baseline, chr >= 2) %>% group_by(tissue) %>% summarise(sum = sum(num)))
  
  
  ### feature
  feature <- feature %>% mutate(tissue = gsub('\\.[0-9]+', '', row.names(.)))
  
  count <- gather(feature, Feature, Count, enhancer_d:synonymous_variant_d) 
  
  # merge baseline and feature
  merge <- merge(count, baseline, by = 'tissue')
  names(merge) <- c('tissue', 'Feature', 'Count', 'Baseline')
  
  return(merge)
}

############# Combine count information
Ctrl.Comb <- function(which.QTL){
  
  
  ############### Observed
  baseline <- read.table(paste0('Result/02_Functional_annotation_enrichment/00_input/ID_RDS/', which.QTL, '.obs.baseline')) %>% mutate(tissue = gsub('\\.[0-9]+', '', row.names(.)))
  feature <- read.table(paste0('Result/02_Functional_annotation_enrichment/00_input/ID_RDS/', which.QTL, '.obs.feature')) %>% mutate(tissue = gsub('\\.[0-9]+', '', row.names(.)))
  rqtl.obs <- Get.obs.Count(baseline, feature)
  
  
  ############## Ctrl
  all <- c()
  for(i in 1:30){
    
    baseline <- read.table(paste0('Result/02_Functional_annotation_enrichment/00_input/ID_RDS/all_', which.QTL, '_random_ID/time.', i, '.baseline')) %>% mutate(tissue = gsub('\\.[0-9]+', '', row.names(.)))
    feature <- read.table(paste0('Result/02_Functional_annotation_enrichment/00_input/ID_RDS/all_', which.QTL, '_random_ID/time.', i, '.feature')) %>% mutate(tissue = gsub('\\.[0-9]+', '', row.names(.)))
    # count and output
    rqtl.ctrl <- Get.Count(baseline, feature)
    write.table(rqtl.ctrl, paste0('Result/02_Functional_annotation_enrichment/00_input/ID_RDS/all_rQTL_random_ID/time.', i, '.count'), row.names = F, quote = F, sep = '\t')
    
    merge <- merge(rqtl.obs, rqtl.ctrl, by = c('tissue', 'Feature'))
    names(merge) <- c('tissue', 'Feature', 'Count.obs', 'Baseline.obs', 'Count.ctrl', 'Baseline.ctrl')
    merge$Time <- i
    
    all <- rbind(all, merge)
  }
  
  return(all)
}

################################################    Main    ##################################################
library(data.table)
library(tidyverse)
setwd('/workspace/rsrch1/ychen/Projects/rQTL/cis_QTL/')


rqtl <- Ctrl.Comb('rQTL')
eqtl <- Ctrl.Comb('eQTL')

rqtl$Label <- 'rhyQTL'
eqtl$Label <- 'eQTL'

# matche tissue names in eQTL with those in rQTL
color <- fread('Data/color_annotation_v3.txt')
for(x in unique(eqtl$tissue)){eqtl$tissue[eqtl$tissue == x] <- color$tissue_rQTL[color$tissue_eQTL == x]}

all <- rbind(rqtl, eqtl)
all$Ratio <- ((all$Count.obs + 0.5)/(all$Count.ctrl + 0.5))/((all$Baseline.obs + 0.5)/(all$Baseline.ctrl + 0.5))
all$Feature <- gsub('_d', '', all$Feature)
all$Feature <- gsub('^X', '', all$Feature)

write.table(all, 'Result/02_Functional_annotation_enrichment/baseline/vep_feature_enrichment.txt', quote = F, sep = '\t', row.names = F)

#######################
#
# generate table S in the paper
#
########################

# names
num_vars <- 30
var_names <- list()
for (i in 1:num_vars) {
  count_name <- paste("Count.ctrl.", i, sep = "")
  baseline_name <- paste("Baseline.ctrl.", i, sep = "")
  ratio_name <- paste("Odds ratio.", i, sep = "")
  
  var_names[[i]] <- c(count_name, baseline_name, ratio_name)
}
var_names <- unlist(var_names)

merge <- c()
for(qtl in unique(all$Label)){
  
  
  combine <- filter(all, Label == qtl & Time == 1)[,1:4]
  combine$QTL.type <- qtl
  combine <- select(combine, tissue, Feature, QTL.type, Count.obs, Baseline.obs)
  tmp <- filter(all, Label == qtl)
  
  
  for(i in 1:30){
    
    combine <- cbind(combine, select(filter(tmp, Time == i), contains(c('ctrl', 'Ratio'))))
    
  }  
  
  names(combine) <- c('Tissue', 'Feature', 'QTL.type', 'Count.obs', 'Baseline.obs', var_names)
  
  merge <- rbind(merge, combine)
}

write.table(merge, 'Result/02_Functional_annotation_enrichment/baseline/vep_feature_enrichment_convert.txt', quote = F, sep = '\t', row.names = F)
convert <- fread('Result/02_Functional_annotation_enrichment/baseline/vep_feature_enrichment_convert.txt') %>% as.data.frame()
names(convert)
for(i in c(as.numeric(grep('Odds', names(convert))))){
  
  convert[,i] <- signif(convert[, i], digit = 2)
  
}
write.table(convert, 'Result/02_Functional_annotation_enrichment/baseline/vep_feature_enrichment_convert_2digit.txt', row.names = F, quote = F, sep = '\t')


############################### PLOT ###############################
library(data.table)
library(tidyverse)
all <- fread('Result/02_Functional_annotation_enrichment/baseline/vep_feature_enrichment.txt')
summary <- as.data.frame(all %>% group_by(tissue, Label, Feature) %>% summarise(median = median(Ratio), sd = sd(Ratio), mean = mean(Ratio)))


theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) 
#theme(legend.position = "none")



median <- as.data.frame(summary %>% group_by(Label, Feature) %>% summarise(median = median(median)))
sd <- as.data.frame(summary %>% group_by(Label, Feature) %>% summarise(sd = sd(median)))
median <- merge(median, sd, by = c('Label', 'Feature'))
median$Label <- factor(median$Label, levels = c('rhyQTL', 'eQTL'))

enrichment.ratio <- as.data.frame(select(median, Feature, Label, median) %>% spread(Label, median) %>% mutate(rhyQTLvseQTL = rhyQTL/eQTL))
enrichment.ratio_sort <- enrichment.ratio[order(-enrichment.ratio$rhyQTLvseQTL), ]
enrichment.ratio_sort$Feature <- factor(enrichment.ratio_sort$Feature, levels = enrichment.ratio_sort$Feature)


p <- ggplot(data = median, aes(x =  Feature, y = median, fill = Label)) + geom_col(position = position_dodge(width = 0.5), width = 0.5) +  theme +
  geom_hline(yintercept = 1, color = 'grey30', linetype = 'dashed') + labs(x = '', y = 'Odds ratio') +
  geom_errorbar(aes(ymin = median - sd, ymax = median + sd), width = 0, position = position_dodge(width = 0.5))+
  labs(y = "Enrichment", x = '') +
  scale_x_discrete(limits = levels(enrichment.ratio_sort$Feature),
                   label = rev(c('Splice donor', 'Promoter', 'Stop gained', 'Frameshift', '5 UTR', 'NC transcript', 'Splice acceptor', '3 UTR',  'Splice region', 'Missense', 'Synonymous',
                             'Intron', 'CTCF binding site', 'TF binding site', 'Promoter-flanking', 'Open chromatin', 'Enhancer'))) +
  scale_fill_manual(values = c('red', '#1C86EE'))  + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + theme(legend.position = "none")
p



p1 <- ggplot(data = enrichment.ratio_sort, aes(x = Feature, y = rhyQTLvseQTL)) + geom_point(size = 0.4) + theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
p1

library(grid)
pdf('Figure/02_Functional_annotation_enrichment/vep_baseline_enrichment.pdf', height = 6, width = 4, useDingbats = F)
print(p, vp=viewport(1, .4, x =.5, y = .7))
print(p1, vp=viewport(1, .5, x =.5, y = .2))
dev.off() 










