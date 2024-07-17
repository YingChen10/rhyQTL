setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')
results <- readRDS('00_rQTL_mapping/04_Result/03_compare_multiple_times_rQTL_pvalue.rds')


results.pval <- mclapply(results, function(tissue.list){
  
  tissue.list$max_amp <- apply(select(tissue.list, contains("amp.c")), 1, function(x) max(x))
  tissue.list$max_amp_idx <- apply(select(tissue.list, contains("amp.c_")), 1, function(x) which.max(x))
  tissue.list$pval <- apply(select(tissue.list, contains("pval_"), max_amp_idx), 1, function(x) (x[x[3]]))
  tissue.list$qval <- apply(select(tissue.list, contains("qval_"), max_amp_idx), 1, function(x) (x[x[3]]))
  
  return(tissue.list)
  
})

###########################################################
# # of rGenes and rQTLs with different p value cutoff
###########################################################

## retrieve all cis-rQTL, use p value cutoff 10-4
pval.cutoff <- c(0.01, 0.001, 0.0001, 0.00001, 0.00001, 0.000001, 0.0000001, 0.00000001)

all <- c()
for(pval.tmp in pval.cutoff){
  
  length.list <- mclapply(results.pval, function(tissue.list){
    
    tissue.list <- filter(tissue.list, pval <= pval.tmp)
    rqtl_length  = length(unique(tissue.list$ID))
    rgene_length  = length(unique(tissue.list$gene))
    
    a <- data.frame(rqtl_length = rqtl_length,
                    rgene_length  = rgene_length)
    return(a)
  })
  
  length <- as.data.frame(do.call(rbind, length.list))
  length$cutoff <- pval.tmp
  
  all <- rbind(all, length)
}

all <- rownames_to_column(all, var = 'tissue')
all$tissue[all$cutoff < 0.01] <- gsub('[0-9]$', "", all$tissue[all$cutoff < 0.01])
all$tissue <- factor(all$tissue, levels = unique(all$tissue))
write.table(all, './00_rQTL_mapping/04_Result/rQTL.rGene.num.p.cutoff.txt', quote = F, sep = '\t', row.names = F)

all <- read.table('./00_rQTL_mapping/04_Result/rQTL.rGene.num.p.cutoff.txt', header = T)
all$tissue <- factor(all$tissue, levels = unique(all$tissue))
color <- fread('Data/color_annotation_v3.txt')


theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) 



p1 <- ggplot(data =  all, aes(x = -log10(cutoff), y = rgene_length, color = tissue)) + geom_line() +  scale_y_log10() +
  labs(y = "# of rGenes", x = '-log10(P value)') +
  scale_color_manual(values = color$ColorPlot[match(levels(all$tissue), color$tissue_rQTL)]) + theme +
  theme(legend.position = "none") + geom_vline(xintercept = 4, color = 'grey', linetype = 'dashed') 
p1

p2 <-  ggplot(data =  all, aes(x = -log10(cutoff), y = rqtl_length, color = tissue)) + geom_line() +  scale_y_log10() +
  labs(y = "# of rQTL", x = '-log10(P value)') +
  scale_color_manual(values = color$ColorPlot[match(levels(all$tissue), color$tissue_rQTL)]) + theme +
  theme(legend.position = "none") + geom_vline(xintercept = 4, color = 'grey', linetype = 'dashed') 
p2  

#########################################################
# # of rgenes and rQTLs with different amplitude cutoff
##########################################################

count <- c()
count <- mclapply(results.pval, function(tissue.list){
  
  tissue.list <- filter(tissue.list, pval <= 0.0001)
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
  
}, mc.cores = 20)

# # of rQTLs and genes involved in all tissues
count_comb <- do.call(rbind, count)
count_comb <- rownames_to_column(count_comb, var = 'tissue')
count_comb$tissue <- gsub('\\.[0-9]+', '', count_comb$tissue)
write.table(count_comb, './00_rQTL_mapping/04_Result/rQTL.rGene.num.Amp.cutoff.txt', quote = F, sep = '\t', row.names = F)

################
# 2.1 plot
############

count_comb <- read.table('./00_rQTL_mapping/04_Result/rQTL.rGene.num.Amp.cutoff.txt', header = T)
count_comb$tissue <- factor(count_comb$tissue, levels = levels(all$tissue))

p3 <- ggplot(filter(count_comb, cutoff >= log2(1.5))) + geom_line(aes(x = cutoff, y = gene_num, colour = tissue)) + scale_y_log10() +
  labs(title = '', y = "# of rGenes", x = "Fold change of peak to trough)") + theme + 
  scale_color_manual(values = color$ColorPlot[match(levels(all$tissue), color$tissue_rQTL)]) +
  theme(legend.position = "none") + geom_vline(xintercept = log2(1.5), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = c(0, log2(1.5), 1, 2, 3, 4, 5, 6), labels = c(1, 1.5, 2, 4, 8, 16, 32, 64))
p3
p4 <- ggplot(filter(count_comb, cutoff >= log2(1.5))) + geom_line(aes(x = cutoff, y = rQTL_num, colour = tissue)) + scale_y_log10() +
  labs(title = '', y = "# of rQTLs", x = "Fold change of peak to trough)") + theme + 
  scale_color_manual(values = color$ColorPlot[match(levels(all$tissue), color$tissue_rQTL)]) +
  theme(legend.position = "none") + geom_vline(xintercept = log2(1.5), linetype = "dashed", color = 'grey') +
  scale_x_continuous(breaks = c(0, log2(1.5), 1, 2, 3, 4, 5, 6), labels = c(1, 1.5, 2, 4, 8, 16, 32, 64))

library(grid)
pdf('Figure/00_rQTL_mapping/rGene.rQTL.number.cutoff.pdf', width = 6, height = 6, useDingbats = F)
print(p1, vp=viewport(.45, .42, x = .7,y = .8))
print(p2, vp=viewport(.45, .42, x = .25,y = .8))
print(p3, vp=viewport(.45, .46, x = .7,y = .3))
print(p4, vp=viewport(.45, .46, x = .25,y = .3))
dev.off()


