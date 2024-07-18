library(data.table)
library(tidyverse)

# combine all h2
file <- list.files('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/LDSC/02_Heritability_top_pval/Liver', pattern = 'results')
all <- c()
for(x in file){
  
  tmp <- fread(paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/LDSC/02_Heritability_top_pval/Liver/', x))
  tmp$trait <- gsub('.ldsc.results', '', x)
  
  all <- rbind(all, tmp)
}
all <- filter(all, Category != 'baseL2_0')
all$Category <- gsub('L2_0', '', all$Category)

write.table(all, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/LDSC/heritability_liver_combine_top_pval.txt', quote = F, sep = '\t', row.names = F)
############################################################################## 
# only keep observed QTLs
all <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/LDSC/heritability_liver_combine_top_pval.txt')
obs <- filter(all, Category == 'rQTL.obs' | Category == 'eQTL.obs')
#### plot observed h2
theme <- theme_bw() +
  theme(axis.line = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        #axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) 

obs$trait.1 <- str_split_fixed(obs$trait, '\\.', 3)[ ,3]

write.table(obs, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/LDSC/observed_QTL_heritability_top_pval.txt', quote = F, sep = '\t', row.names = F)
obs <-  fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/LDSC/observed_QTL_heritability_top_pval.txt')
obs$Category <- gsub('.obs', '', obs$Category)
obs$Category <- factor(obs$Category, levels = rev(unique(obs$Category)))
obs$trait.1 <- gsub('-in-chylomicrons-and-extremely-large-VLDL-UKB-data-field-23485', '', obs$trait.1) %>% gsub('-levels', '', .) %>%
  gsub('Concentration-of-', '', .) %>% gsub('Average-diameter-for-', "", .) %>% gsub('-UKB-data-field-23405', '', .)
obs$trait.1[obs$trait.1 == 'Ratio-of-saturated-fatty-acids-to-total-fatty-acids'] <- 'Saturated vs total fatty acids'
obs$trait.1 <- gsub('-', ' ', obs$trait.1)
obs$trait.1 <- factor(obs$trait.1, levels = rev(sort(unique(obs$trait.1))))

# Proportion of heitablity
p1 <- ggplot(data = obs, aes(x = trait.1, y = Prop._h2, fill = Category)) + geom_col(position = "dodge", width = 0.5) + 
  coord_flip() + theme + ylab('Proportion of heitablity') +
  scale_fill_manual(values = c('#1C86EE', 'red')) +
  geom_errorbar(aes(ymin = Prop._h2 - Prop._h2_std_error * 1.96, ymax = Prop._h2 + Prop._h2_std_error * 1.96), width = 0,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = median(obs$Prop._h2[obs$Category == 'eQTL']), color = '#1C86EE') + 
  geom_hline(yintercept = median(obs$Prop._h2[obs$Category == 'rQTL']), color = 'red') +
  theme(legend.position = "none")

p1

# tau
p2 <- ggplot(data = obs, aes(x = trait.1, y = Coefficient/1e-06, fill = Category)) + geom_col(position = "dodge", width = 0.5) + 
  coord_flip() + theme +  ylab(bquote(italic(tau)^{"*"} ~ "(1e-06)")) +
  scale_fill_manual(values = c('#1C86EE', 'red')) +
  geom_errorbar(aes(ymin = (Coefficient-1.96 * Coefficient_std_error)/1e-06, ymax = (Coefficient+1.96 * Coefficient_std_error)/1e-06), 
                width = 0,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = median(obs$Coefficient[obs$Category == 'eQTL'])/1e-06, color = '#1C86EE') + 
  geom_hline(yintercept = median(obs$Coefficient[obs$Category == 'rQTL'])/1e-06, color = 'red') +
  theme(legend.position = "none")

p2

# enrichment
p3 <- ggplot(data = obs, aes(x = trait.1, y = Enrichment, fill = Category)) + geom_col(position = "dodge", width = 0.5) + 
  coord_flip() + theme +  ylab(bquote(paste("Pr(",italic(h)[g]^{2},")/Pr(SNPs)",sep=""))) +
  scale_fill_manual(values = c('#1C86EE', 'red')) +
  geom_errorbar(aes(ymin = Enrichment-1.96 * Enrichment_std_error, ymax = Enrichment+1.96 * Enrichment_std_error), 
                width = 0,
                position = position_dodge(width = 0.5)) +
  geom_hline(yintercept = median(obs$Enrichment[obs$Category == 'eQTL']), color = '#1C86EE') + 
  geom_hline(yintercept = median(obs$Enrichment[obs$Category == 'rQTL']), color = 'red') +
  theme(legend.position = "none")

p3


pdf('Figure/03_enrich_GWAS/lsdc_top_pval.pdf', width = 18, height = 2.5, useDingbats = F)
print(p1, vp=viewport(.15, .95, x = .2, y = .5))
print(p2, vp=viewport(.15, .95, x = .45, y = .5))
print(p3, vp=viewport(.15, .95, x = .7, y = .5))
dev.off()

median(obs$Prop._h2[obs$Category == 'eQTL'])
median(obs$Prop._h2[obs$Category == 'rQTL'])




