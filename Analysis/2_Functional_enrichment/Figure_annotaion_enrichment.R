
################################################    Main    ##################################################
library(data.table)
library(tidyverse)
setwd('/workspace/rsrch1/ychen/Projects/rQTL/cis_QTL/')



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










