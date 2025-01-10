setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL')

tissue.list <- list.files('13_Results/00_all_rhyQTL/') %>% gsub('.txt', '', .)

###########################################################
# # of rGenes and rQTLs with different p value cutoff
###########################################################

## retrieve all cis-rQTL, use p value cutoff 10-4
pval.cutoff <- c(0.01, 0.001, 0.0005, 0.0001, 0.00001, 0.00001, 0.000001, 0.0000001, 0.00000001)

all <- c()
for(x in tissue.list){
  
  message(x)
  
  # observed rhyQTLs
  tmp <- readRDS(paste0(dir, cat, '00_rQTL_mapping/13_combine/', x, '.rds')) 
  
  for(pval.tmp in pval.cutoff){
    
    tissue.list <- tmp %>%  filter(., pval < pval.tmp)
    tissue.list$p_adj <- p.adjust(tissue.list$HANOVA.norm, method = "BH")
    tissue.list <- filter(tissue.list, p_adj < 0.01)
    tmp.p.cutoff <- unique(select(tissue.list, gene, ID))
    
    all <- data.frame(
      
      tissue = x, 
      rgene_length = length(unique(tmp.p.cutoff$gene)),
      rqtl_length = length(unique(tmp.p.cutoff$ID)),
      cutoff = pval.tmp
      
    ) %>% rbind(all, .)
    
  }
}  

all$tissue <- factor(all$tissue, levels = unique(all$tissue))
write.table(all, '13_Results/rhyQTL.rhyGene.num.p.cutoff.txt', quote = F, sep = '\t', row.names = F)

all <- read.table('13_Results/rhyQTL.rhyGene.num.p.cutoff.txt', header = T)
all$tissue <- factor(all$tissue, levels = unique(all$tissue))
color <- fread('Data/color_annotation_v3.txt')


theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) 



p1 <- ggplot(data =  all, aes(x = -log10(cutoff), y = rgene_length, color = tissue)) + geom_line() +  scale_y_log10() +
  labs(y = "# of rhyGenes", x = '-log10(P value)') +
  scale_color_manual(values = color$ColorPlot[match(levels(all$tissue), color$tissue_rQTL)]) + theme +
  theme(legend.position = "none") + geom_vline(xintercept = -log10(5e-4), color = 'grey', linetype = 'dashed') 
p1

p2 <-  ggplot(data =  all, aes(x = -log10(cutoff), y = rqtl_length, color = tissue)) + geom_line() +  scale_y_log10() +
  labs(y = "# of rhyQTL", x = '-log10(P value)') +
  scale_color_manual(values = color$ColorPlot[match(levels(all$tissue), color$tissue_rQTL)]) + theme +
  theme(legend.position = "none") + geom_vline(xintercept = -log10(5e-4), color = 'grey', linetype = 'dashed') +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000))
p2  

#########################################################
# # of rgenes and rQTLs with different amplitude cutoff
##########################################################
tissue.list <- list.files('13_Results/00_all_rhyQTL/') %>% gsub('.txt', '', .)

results.pval <- lapply(tissue.list, function(tissue){
  
  tmp <- fread(paste0('13_Results/00_all_rhyQTL/', tissue, '.txt'))
  return(tmp)
  
})
names(results.pval) <- tissue.list

count <- c()
count <- lapply(results.pval, function(tissue.list){
  
  tissue.list <- filter(tissue.list, pval <= 5e-4)
  # set break
  breaks <- seq(0, max(max(tissue.list$max_amp),1), by=0.1)
  breaks <- c(breaks, log2(1.5))
  
  # subset rQTLs under each break cutoff
  subset_list <- lapply(breaks, function(b) {
    tissue.list[tissue.list$max_amp > b, c('ID', 'gene')]
  })
  
  # get the # of rQTLs and genes
  count_breaks <- lapply(subset_list, function(sub_results) {
    rQTL_num <- length(unique(sub_results$ID))
    gene_num <- length(unique(sub_results$gene))
    return(data.frame(rQTL_num = rQTL_num, gene_num = gene_num))
  })
  count_break <- do.call(rbind, count_breaks)
  count_break$cutoff <- breaks
  return(count_break)
  
})

# # of rQTLs and genes involved in all tissues
count_comb <- do.call(rbind, count)
count_comb <- rownames_to_column(count_comb, var = 'tissue')
count_comb$tissue <- gsub('\\.[0-9]+', '', count_comb$tissue)
write.table(count_comb, '13_Results/rhyQTL.rhyGene.num.Amp.cutoff.txt', quote = F, sep = '\t', row.names = F)

################
# 2.1 plot
############
count_comb <- read.table('13_Results/rhyQTL.rhyGene.num.Amp.cutoff.txt', header = T)
count_comb$tissue <- factor(count_comb$tissue, levels = unique(count_comb$tissue))

p3 <- ggplot(filter(count_comb, cutoff >= log2(1.5))) + geom_line(aes(x = cutoff, y = gene_num, colour = tissue)) + 
  scale_y_log10() +
  labs(title = '', y = "# of rhyGenes", x = "Fold change of peak to trough)") + theme + 
  scale_color_manual(values = color$ColorPlot[match(levels(count_comb$tissue), color$tissue_rQTL)]) +
  theme(legend.position = "none") + geom_vline(xintercept = log2(1.5), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = c(0, log2(1.5), 1, 2, 3, 4, 5, 6), labels = c(1, 1.5, 2, 4, 8, 16, 32, 64))
p3
p4 <- ggplot(filter(count_comb, cutoff >= log2(1.5))) + geom_line(aes(x = cutoff, y = rQTL_num, colour = tissue)) + scale_y_log10() +
  labs(title = '', y = "# of rhyQTLs", x = "Fold change of peak to trough)") + theme + 
  scale_color_manual(values = color$ColorPlot[match(levels(count_comb$tissue), color$tissue_rQTL)]) +
  theme(legend.position = "none") + geom_vline(xintercept = log2(1.5), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = c(0, log2(1.5), 1, 2, 3, 4, 5, 6), labels = c(1, 1.5, 2, 4, 8, 16, 32, 64)) +
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000))
p4
library(grid)
pdf('13_Figures/00_rhyGene.rhyQTL.number.cutoff.pdf', width = 6, height = 6, useDingbats = F)
print(p1, vp=viewport(.45, .415, x = .7,y = .7))
print(p2, vp=viewport(.45, .415, x = .25,y = .7))
print(p3, vp=viewport(.45, .46, x = .7,y = .3))
print(p4, vp=viewport(.45, .46, x = .25,y = .3))
dev.off()

