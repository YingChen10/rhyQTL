



### get the number of tissues in  
# cd /workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/00_rQTL_mapping/00_Meff_input
# awk '{key = $1 "\t" $2; count[key]++} END {for (k in count) print count[k], k}' *.txt > ../../13_Results/02_tissue_sharing_rhyGene/all_tissue_gene_snp_pair.txt


library(data.table)
library(tidyverse)
dir <- '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/'


file_path <- paste0(dir, '02_tissue_sharing_rhyGene/all_tissue_gene_snp_pair.txt')
all <- read.table(file_path, header = FALSE, sep = "", stringsAsFactors = FALSE)
names(all) <- c('num', 'id', 'gene')
all$idx <- paste0(all$gene, ':', all$id)


tissue.specific <- fread(paste0(dir, '02_tissue_sharing_rhyGene/tissue_specific_rhyGene.txt'))
head(tissue.specific)

combine <- data.frame()
for(x in unique(tissue.specific$tissue)){
  
  tmp <- fread(paste0(dir,'00_all_rhyQTL/', x, '.txt'))
  tmp.filter <- filter(tmp, gene %in% tissue.specific$rhyGene[tissue.specific$tissue == x])
  
  index <- paste0(tmp.filter$gene, ':', tmp.filter$ID)
  
  combine <- data.frame(tissue = x, idx = index) %>% rbind(combine, .)
  
}
write.table(combine, paste0(dir, '02_tissue_sharing_rhyGene/tissue_specific_rhyGene_rhyQTL_pair.txt'), quote = F, sep = '\t', row.names = F)
head(all)
head(combine)

filter <- merge(all, combine, by = 'idx')
write.table(filter, paste0(dir, '02_tissue_sharing_rhyGene/tissue_specific_rhyGene_rhyQTL_pair_num.txt'), quote = F, sep = '\t', row.names = F)



filter <- fread(paste0(dir, '02_tissue_sharing_rhyGene/tissue_specific_rhyGene_rhyQTL_pair_num.txt'))
#median <- as.data.frame(filter %>% group_by(tissue) %>% summarise(median = median(num)))


result <- filter %>%
  group_by(tissue, gene) %>%
  filter(num == max(num)) %>%
  slice_sample(n = 1) %>%  
  ungroup()

median <- as.data.frame(result %>% group_by(tissue) %>% summarise(median = median(num)))


theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) 


color <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Data/color_annotation_v3.txt')
sample.size <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Data/00_sample_size_with_phenotype.txt')
median <- merge(median, color, by.x = 'tissue', by.y = 'tissue_rQTL')
median <- merge(median, sample.size, by = 'tissue')
median <- median[order(-median$median),]
median$ShortName <- factor(median$ShortName, levels = median$ShortName)


head(median)
str(median)

p1 <- ggplot(data = median, aes(x = ShortName, y = median)) + 
  geom_col(color = 'black', fill = 'skyblue') + theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
  labs(x = '', y = 'Median of tissue numbers tested') 
p1

# heatmap of sample size
p2 <- ggplot(data = median, aes(x = ShortName, y = 1, fill = samplesize)) + geom_tile() + 
  scale_fill_gradientn(colours = c("yellow", "orange", "darkred")) + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))

p2 

p3 <- ggplot(data = median, aes(x = ShortName, y = 1, fill = ShortName)) + geom_tile() + 
  scale_fill_manual(values = median$ColorPlot[match(levels(median$ShortName), median$ShortName)]) + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) + 
  theme(legend.position = "none")
p3 

library(grid)
pdf('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Figures/tissue_specific.pdf', height = 9, width = 4, useDingbats = F)
print(p1, vp=viewport(1, .3, x = .5, y = .8))
print(p2, vp=viewport(1, .3, x = .5, y = .5))
print(p3, vp=viewport(1, .3, x = .5, y = .2))
dev.off()  


