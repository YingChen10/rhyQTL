library(data.table)
library(parallel)
library(tidyverse)

setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')
color <- fread('Data/color_annotation_v3.txt')

tissue.list <- list.files('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/00_all_rhyQTL') %>% gsub('.txt', '', .)



#library(ggvenn)
#pdf("./Figures/01_eGene_sGene_overlap/eQTL_overlap/rqtl_eqtl_overlap.pdf", width = 3.6, height = 3.6, useDingbats = F)
overlap.list <- lapply(tissue.list, function(x){
  
  rqtl <- fread(paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/00_all_rhyQTL/', x, '.txt'))
  eqtl_tissue <- color$tissue_eQTL[color$tissue_rQTL == x]
  tmp <- fread(paste0('../eQTLs/GTEx_Analysis_v8_eQTL/', eqtl_tissue, '.v8.signif_variant_gene_pairs.txt.gz'))
  
  rqtl_id <- unique(rqtl$ID)
  eqtl_id <- unique(tmp$variant_id)
  
  # y <- list(rQTL = rqtl_id, eQTL = eqtl_id)
  
  # print(ggvenn(y, fill_color = c("#0073C2FF", "#EFC000FF"),stroke_size = 0.5, set_name_size = 4, show_percentage = FALSE),vp = viewport(1, 1, x=.5, y=.5))
  # grid.text(x, x = unit(0.5, "npc"), y = unit(0.9, "npc"), gp = gpar(fontsize = 12, fontface = "bold"))
  # 
  # grid.newpage()
  
  overlap_result <- data.frame(
    tissue = x,
    common = length(intersect(rqtl_id,  eqtl_id)),
    eqtl_count = length(unique(eqtl_id)),
    rqtl_count = length(unique(rqtl_id)))
  
  return(overlap_result)
  
})
#dev.off()


overlap <- as.data.frame(do.call(rbind, overlap.list))
overlap$rqtl_common_frac <- overlap$common/overlap$rqtl_count
overlap$rqtl_uniq_frac <- 1 - overlap$rqtl_common_frac

color <- fread('./Data/color_annotation_v3.txt')
overlap <- merge(color, overlap, by.x = 'tissue_rQTL', by.y = 'tissue')
overlap <- overlap[order(overlap$rqtl_common_frac),]
overlap$ShortName <- factor(overlap$ShortName, levels = overlap$ShortName)
write.table(overlap, './13_Results/02_rQTL_eQTL_overlap.txt', quote = F, sep = '\t', row.names = F)

###########
# plot QTL overlap
###########
overlap <- fread('./13_Results/02_rQTL_eQTL_overlap.txt')
overlap <- overlap[order(overlap$rqtl_common_frac),]
head(overlap)
range(overlap$rqtl_common_frac)

overlap$ShortName <- factor(overlap$ShortName, levels = overlap$ShortName)
head(overlap)
overlap_convert <- gather(select(overlap, ShortName, rqtl_common_frac, rqtl_uniq_frac), type, fraction, rqtl_common_frac, rqtl_uniq_frac)
overlap_convert$ShortName <- factor(overlap_convert$ShortName, levels = overlap$ShortName)

theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.ticks.x = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        #axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) + 
  theme(legend.position = "none")

p1 <- ggplot(data = overlap_convert, aes(x = ShortName, y = fraction, fill = type)) + geom_col() +
  scale_fill_manual(values = c('#1C86EE', 'red')) + theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + ylab('Fraction of rQTLs')
p1


p2 <- ggplot(data = overlap, aes(x = ShortName, y = 1, fill = ShortName)) + geom_col() + theme + 
  scale_fill_manual(values = overlap$ColorPlot) +
  theme(legend.position = "none") + labs(x = '', y = '')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
p2

library(grid)
pdf('./13_Figures/02_eQTL_rQTL_overlap.pdf', height = 4, width = 4, useDingbats = F)
grid.arrange(p1, p2, ncol = 1)
dev.off()  


