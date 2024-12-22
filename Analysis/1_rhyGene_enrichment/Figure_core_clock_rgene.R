library(data.table)
library('tidyverse')
setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')
core.clock <- c('ARNTL', 'RORC', 'CRY1', 'CRY2', 'NPAS2', 'NR1D1', 'NR1D2', 'PER1', 'PER2', 'PER3', 'CLOCK')

rgene <- fread('./13_Results/00_tissue_rhyGene.name.list')
names(rgene) <- c('rhy.Gene', 'name', 'tissue')
annot <- fread('all_gene.annot')

tissue <- unique(rgene$tissue)
all <- c()

# annotate weather the core clock gene is rhyGene or not 
for(x in core.clock){
  
  tmp <- annot$V5[annot$V7 == x]
  result <- data.frame(tissue = tissue, rGene = tmp, name = x, hl = 0)
  result$hl[match(filter(rgene, rhy.Gene == tmp)$tissue, result$tissue)]  <- 1
  all <- rbind(all, result)
  
}

theme <- theme_bw() +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.title = element_text(size = 8)) + 
  theme(legend.position = "none")


color <- fread('Data/color_annotation_v3.txt')

all$tissue <- factor(all$tissue, levels = unique(all$tissue))

p1 <- ggplot(all, aes(x = tissue, y = name, fill = as.character(hl))) + 
  geom_tile(colour = "black") + 
  theme + 
  scale_fill_manual(values = c("grey100", "red")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  labs(x = '', y = "") + scale_y_discrete(limits  = rev(sort(unique(all$name))))
p1
# 
# p1 <- ggplot(all, aes(x = tissue, y = name, fill = as.character(hl))) + 
#   geom_tile(colour = "black") + 
#   theme + 
#   scale_fill_manual(values = c("grey100", "orange")) +
#   theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
#   labs(x = '', y = "") + scale_y_discrete(limits  = rev(sort(unique(all$name))))
# p1

color <- filter(color, tissue_rQTL %in% levels(all$tissue))
p2 <- ggplot(color, aes(x = ShortName, y = 1, fill = tissue_rQTL)) + 
  geom_tile(colour = "white") + 
  theme + 
  scale_fill_manual(values = color$ColorPlot[match(levels(all$tissue), color$tissue_rQTL)])+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) +
  labs(x = '', y = "") 
p2

library(grid)
pdf('./13_Figures/Figure_core_clock_rgene.pdf', width = 6, height = 9, useDingbats = T)
print(p1, vp=viewport(1, .42, x = .5, y = .7))
print(p2, vp=viewport(1, .35, x = .5, y = .2))
dev.off()

