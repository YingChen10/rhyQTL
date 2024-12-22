tissue.file <- list.files('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/00_all_rhyQTL/', pattern = 'txt')
tissue.list <- gsub('.txt', '', tissue.file)
#tissue.list <- tissue.list[tissue.list != 'Kidney-Cortex']
file.list <- list.files(paste0(dir, cat, 'split_pos'))

all <- data.frame()
for(tissue in tissue.list){
  
  
  message(paste0(tissue, ' Start!'))
  
  result <- lapply(file.list, function(file){
    
    file.tmp <- fread(paste0(dir, cat, '00_rQTL_mapping/02_Rhythm_compare/', tissue, '/', file))
    #compare <-  fread(paste0(dir, cat, '00_rQTL_mapping/03_Rhythm_compare_multiple_times/', tissue, '/', file, '.filter'))
    return(file.tmp)
})
  result <- as.data.frame(do.call(rbind, result))
  names(result)[names(result) == 'Gene'] <- 'gene'
  rqtl <- fread(paste0(dir, cat, '13_Results/00_all_rhyQTL/', tissue, '.txt'))  
  
  result <- merge(result, rqtl, by = c('gene', 'ID')) %>%
    select(chosen_model.x, chosen_model_mean) 
  
  a <- as.data.frame(table(result$chosen_model_mean)) %>% mutate(tissue = tissue)
  
  all <- rbind(all, a)
  
}


all.convert <- spread(all, Var1, Freq)

names(all.convert) <- c('tissue','same', 'diff')
all.convert$ratio <- all.convert$diff/(all.convert$same + all.convert$diff)
all.convert
write.table(all.convert, '13_Results/rhyQTL_diff_expr.txt', quote = F, sep = '\t', row.names = F)

all.convert <- fread('13_Results/rhyQTL_diff_expr.txt')
count <- fread('13_Results/00_rhyGene.rhyQTL.num.txt')

# get tissue order in the plot
all.convert <- merge(all.convert, count, by = 'tissue')
all.convert <- all.convert[order(all.convert$ratio),]
all.convert$ShortName <- factor(all.convert$ShortName, levels = all.convert$ShortName)
all.convert$tissue <- factor(all.convert$tissue, levels = all.convert$tissue)
all.convert$ColorPlot <- factor(all.convert$ColorPlot, levels = all.convert$ColorPlot)

theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) +
  theme(legend.position = "none") 

head(all.convert)

# plot: extend each tissue
names(all.convert)[names(all.convert) == 'ratio'] <- 'common'
all.convert$unique <- 1 - all.convert$common
overlap <- gather(all.convert, type, fraction, common, unique)



p1 <- ggplot(data = overlap, aes(x = ShortName, y = fraction, fill = type)) + geom_col() +
  scale_fill_manual(values = c('#1C86EE', 'red')) + theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)) + ylab('Fraction of rQTLs')
p1


p2 <- ggplot(data = overlap, aes(x = ShortName, y = 1, fill = ShortName)) + geom_col() + theme + 
  scale_fill_manual(values = levels(overlap$ColorPlot)) +
  theme(legend.position = "none") + labs(x = '', y = '')+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1))
p2

library(gridExtra)
library(grid)
pdf('./13_Figures/rhyQTL_diff_expr.pdf', height = 4, width = 4, useDingbats = F)
grid.arrange(p1, p2, ncol = 1)
dev.off()  

range(overlap$fraction[overlap$type == 'common'])

overlap.eqtl <- fread('13_Results/02_rQTL_eQTL_overlap.txt')
all.convert <- fread('13_Results/rhyQTL_diff_expr.txt')
head(overlap.eqtl)
head(all.convert)

all <- merge(all.convert, overlap.eqtl, by.x = 'tissue', by.y = 'tissue_rQTL')
ggplot(data = all, aes(x =  ratio, y = rqtl_common_frac)) + geom_point()
cor.test(all$ratio, all$common)

### combined boxplot
p1 <- ggplot() + 
  geom_boxplot(data = all.convert, aes(x = as.character(1), y = ratio * 100), outliers = F) +
  geom_jitter(data = all.convert, aes(x = as.character(1), y = ratio * 100, color = tissue), width = 0.2) +
  theme +
  scale_colour_manual(values = levels(all.convert$ColorPlot)) + 
  labs(x = '', y = '% of rhyQTLs that also affect the mean expression levels between genotypes')
p1

ggplot() + 
  geom_col(data = all.convert, aes(y = tissue, x = ratio * 100, fill = tissue)) +
  theme +
  scale_fill_manual(values = levels(all.convert$ColorPlot)) + 
  labs(x = '', y = '% of rhyQTLs that also affect the mean expression levels between genotypes')

pdf('13_Figures/08_eQTL_overlap/diff_expr_ratio.pdf', width = 2.5, height = 3)
print(p1)
dev.off()


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
overlap_convert <- overlap_convert[order(overlap_convert$ratio),]

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

