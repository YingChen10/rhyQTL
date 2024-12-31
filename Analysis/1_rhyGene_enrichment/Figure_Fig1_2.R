rm(list=ls())
library('tidyverse')
library('data.table')
library('parallel')

setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/00_rQTL_mapping/')
tissue.list <- list.dirs('00_Genotype')[-1] %>% gsub('00_Genotype/', '', .)

###############
#
# get numbers in Figure 1
#
###############
library(paralle)
file <- list.dirs(paste0(dir, cat, '00_rQTL_mapping/00_Genotype'))
tissue <- gsub(paste0(dir, cat, '00_rQTL_mapping/00_Genotype/'), '', file)[-1]


file.list <- list.files(paste0(dir, cat, 'split_pos'))

all <- mclapply(tissue, function(tissue.tmp){
  
  count <- lapply(file.list, function(file){
    
    tmp <- readRDS(paste0(dir, cat, '00_rQTL_mapping/01_Rhythm_regression/', tissue.tmp, '/', file, '.rds'))
    tmp <- as.data.frame(do.call(rbind, tmp)) %>% 
      rownames_to_column(., var = 'idx') %>% 
      mutate(ID = paste0('chr', str_split_fixed(idx, '.chr', 2)[ ,2])) %>%
      mutate(gene = str_split_fixed(idx, '.chr', 2)[ ,1]) 
    
    tmp$max_amp <- apply(select(tmp, contains("amp.")), 1, function(x) max(x))
    tmp$max_amp_idx <- apply(select(tmp, contains("amp.")), 1, function(x) which.max(x))
    tmp$pval <- apply(select(tmp, contains("pval_"), max_amp_idx), 1, function(x) (x[x[3]]))
    filter <- filter(tmp,  pval < 5e-4)
    
    return(
      data.frame(file = file,
                 gene.num = length(unique(filter$gene)),
                 id.num = length(unique(filter$ID)),
                 pair.num = length(unique(filter$idx)))
    )
    
  })
  count <- as.data.frame(do.call(rbind, count))
  count$tissue <- tissue.tmp
  
  return(count)
  
}, mc.cores = 5)

all <- as.data.frame(do.call(rbind, all))

combine <- as.data.frame(all %>% group_by(tissue) %>% summarise(gene.num = sum(gene.num), id.num = sum(id.num), pair.num = sum(pair.num)))
head(combine)

median(combine$gene.num)
# 3231
median(combine$pair.num)
# 620356

#####################################################################################
#
# expressed gene in each tissue as background of enrichment analysis
#
#####################################################################################
file <- list.files('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/GTEx_nor_expression')
file <- file[-match("Kidney-Cortex.txt", file)]

outdir <- paste0(dir, cat, '13_Results/01_expr_gene_bck/') 
if (!file.exists(outdir)) {dir.create(outdir, recursive = TRUE)}

for(x in file){
  
  tmp <- fread(paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/GTEx_nor_expression/', x))
  write.table(unique(str_split_fixed(tmp$EnsemblID, '_', 2)[ ,2]),
              paste0(outdir, gsub('.txt', '', x), '.bck.txt'),
              quote = F, sep = '\t', row.names = F, col.names = F)
  
}


######################################################################################################
#
#                           Figure 1. rhyGenes overlap with identified previously rhythmic genes 
#
#####################################################################################################

###########################################################
# 1. calculate previously identified circadian gene number
###########################################################

file <- list.files('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/GTEx_nor_expression/', pattern = 'txt')

tissue <- gsub('.txt', '', file)
tissue <- tissue[tissue != 'Kidney-Cortex']

dir <- '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/'
cir_gene <- c()
sample_size <- c()
for(x in tissue){
  # input expression file
  expression <- fread(paste0(dir,'GTEx_nor_expression/', x, '.txt'))
  expression <- as.data.frame(expression)
  sample_size <- data.frame(tissue = x, size = length(names(expression[,2:(ncol(expression)-34)]))) %>% rbind(sample_size,.)
  expression <- expression[,c(1, (ncol(expression) - 33) : ncol(expression))]
  expression$EnsemblID <- sapply(strsplit(expression$EnsemblID, '_'), "[", 1)
  expression <- filter(expression, qval.ALL < 0.2 & amp.ALL > 0.5)
  cir_gene <- data.frame(gene = expression$EnsemblID, tissue = x) %>% rbind(cir_gene, .)
}

cir_gene_num <- as.data.frame(table(cir_gene$tissue))
names(cir_gene_num) <- c('tissue', 'cir_gene_number')
merge <- merge(cir_gene_num, sample_size, by = 'tissue')
write.table(merge, '01_eGene_sGene_overlap/rGene_circadian_gene_overlap/circadian_gene_number_from_Science.txt', quote = F, sep = '\t', row.names = F)
write.table(cir_gene, '01_eGene_sGene_overlap/rGene_circadian_gene_overlap/circadian_gene_list.txt', quote = F, sep = '\t', row.names = F)

setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')

cir_gene <- read.table('01_eGene_sGene_overlap/rGene_circadian_gene_overlap/circadian_gene_list.txt', header = T)
annot <- read.table('../gencode.v26.GRCh38.genes.annot')
annot <- unique(annot[, c('V6', "V7")])
names(annot) <- c('gene', 'name')
cir_gene <- merge(cir_gene, annot, by = 'gene')
write.table(cir_gene, '01_eGene_sGene_overlap/rGene_circadian_gene_overlap/circadian_gene_names.list.txt', sep = '\t', quote = F, row.names = F)

#######################################################
#
# output overlap gene and perform enrichment analysis
#
#######################################################
setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')
rgene <- read.table('13_Results/00_tissue_rhyGene.name.list', header = T)
cir_gene <- read.table('01_eGene_sGene_overlap/rGene_circadian_gene_overlap/circadian_gene_names.list.txt', header = T)

outdir <- paste0(dir, cat, '13_Results/01_rhyGene_circadian_gene_overlap/')
if (!file.exists(outdir)) {dir.create(outdir, recursive = TRUE)}

count <- c()
for (x in unique(cir_gene$tissue)){
  
  r.gene <- rgene$name[rgene$tissue==x]
  c.gene <- cir_gene$name[cir_gene$tissue == x]
  
  common <- intersect(r.gene, c.gene)
  r.unique <- setdiff(r.gene, common)
  c.unique <- setdiff(c.gene, common)
  
  outdir.tmp <- paste0(outdir, x ,'/')
  dir.create(outdir.tmp, recursive = TRUE)
  
  write.table(common, paste0(outdir.tmp, 'common.txt'), quote = F, row.names = F, col.names = F)
  write.table(r.unique, paste0(outdir.tmp, 'rhyGene.uniq.txt'), quote = F, row.names = F, col.names = F)
  write.table(c.unique, paste0(outdir.tmp, 'circaGene.uniq.txt'), quote = F, row.names = F, col.names = F)
  
  data <- data.frame(common = paste(strsplit(common, " "), collapse = ","), 
                     r.uniq = paste(strsplit(r.unique, " "), collapse = ","),
                     c.uniq = paste(strsplit(c.gene, " "), collapse = ","))
  
  data <- t(data)
  write.table(data, paste0(outdir.tmp, 'metascape.input.txt'), quote = F, col.names = F, sep = '\t')
  
  
  # get the number of different types of genes
  count <- data.frame(
    tissue = x,
    rhy.gene.specific = length(c.unique),
    rgene.specific = length(r.unique),
    common = length(common)
  ) %>% rbind(., count)
  
  
}


# calculate percentage
count$percent.new.rhy.gene <- count$rgene.specific/(count$rhy.gene.specific+count$rgene.specific+count$common) * 100
count$percent.rhy.gene.genotype.specific <- count$rgene.specific/(count$rhy.gene.specific+count$common) 
median(count$percent.rhy.gene.genotype.specific)

write.table(count, paste0(outdir, '01_rhyGene.cirGene.count.txt'), quote = F, sep = '\t', row.names = F)

library(data.table)
count <- fread(paste0(outdir, '01_rhyGene.cirGene.count.txt'))
## 3 fold more new rhythmic genes were detected in specific genotypes 
count$new_rhythmic_genotype <- count$rgene.specific/(count$rhy.gene.specific + count$common)

median(count$new_rhythmic_genotype)



count <- count[order(-count$rgene.specific),]
count_convert <- gather(count[,1:4], type, number, rhy.gene.specific:common)
head(count_convert)

theme <-  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2)) + 
  theme(legend.position = "none")


count_convert$tissue <- factor(count_convert$tissue, levels =  count$tissue)
count_convert$type <- factor(count_convert$type, levels = c("rgene.specific", "common", "rhy.gene.specific"))
pp1 <- ggplot(data = count_convert, aes(x = tissue, y = number, fill = type)) + geom_col(position = "dodge", width = 0.8) + 
  labs(y = "# of rGenes", x = "# of circadian genes") + theme +
  scale_fill_manual(values = c('#68CAD7', '#F7941D', '#F284AD')) + 
  theme(legend.position = "none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
pp1


color <- fread('Data/color_annotation_v3.txt')
color$tissue_rQTL <- factor(color$tissue_rQTL, levels = levels(count_convert$tissue))
color$FullName <- factor(color$FullName, levels = color$FullName[match(levels(color$tissue_rQTL), color$tissue_rQTL)])

p <- ggplot(data = color, aes(x = FullName, y = 1, fill = FullName)) + geom_tile() +
  scale_fill_manual(values = color$ColorPlot[match(levels(color$FullName), color$FullName)]) + theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
p

plot.outdir <- paste0(dir, cat, '13_Figures/01_rhyGene_circadian_gene_overlap/')
if (!file.exists(plot.outdir)) {dir.create(plot.outdir, recursive = TRUE)}
library('grid')
pdf(paste0(plot.outdir,'rhyGene_rhyGene_overlap.pdf'), width = 5, height = 3.5)
print(p)
print(pp1)
dev.off()


######################################################################################################
#
#                         Figure 2a. # of rhyQTLs and rhyGenes
#
#####################################################################################################
#rm(list = ls()) 
library('tidyverse')
library(parallel)
library('data.table')
setwd('~/Projects/Project03_human_circadian/rQTL/cis_QTL/')

dir <- '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/'
cat <- 'cis_QTL/'

file <- list.dirs(paste0(dir, cat, '00_rQTL_mapping/00_Genotype'))
tissue <- gsub(paste0(dir, cat, '00_rQTL_mapping/00_Genotype/'), '', file)[-1]

num <- c()
for(x in tissue){
  
  tmp <- fread(paste0(dir, cat, '13_Results/00_all_rhyQTL/', x, '.txt'))
  num <- data.frame(
    tissue = x,
    rhyQTL.num = length(unique(tmp$ID)),
    rhyGene.num = length(unique(tmp$gene)),
    pair.num = nrow(tmp)) %>% rbind(num, .)
  
}

# get numbers in Fig1 
median(num$rhyQTL.num)
# 51364
median(num$rhyGene.num)
# 2047
median(num$pair.num)
# 52485

color <- fread('Data/color_annotation_v3.txt')
count <- merge(num, color, by.x = 'tissue', by.y = 'tissue_rQTL')
count <- count[order(count$rhyGene.num),]
count$tissue <- factor(count$tissue, levels = count$tissue)
count$ShortName <- factor(count$ShortName, levels = count$ShortName)

write.table(count, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/00_rhyGene.rhyQTL.num.txt', quote = F, sep = '\t', row.names = F)


library(ggrepel)
library(patchwork)

theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8))

####################################### number of rGene
p1 <- ggplot(data = count, aes(x = ShortName, y = rhyGene.num, fill = tissue)) + geom_col(width = 0.5) + theme_bw() +
  coord_flip() + 
  theme(axis.line = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) + 
  labs(y = "", x = 'Tissues in GTEx') +
  scale_fill_manual(values = count$ColorPlot) +
  theme(legend.position = "none")
p1


############## # of rQTL
p2 <- ggplot(data = count, aes(x = ShortName, y = rhyQTL.num, fill = tissue)) + geom_col(width = 0.5) + theme_bw() +
  coord_flip() + 
  theme(axis.line = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) + 
  labs(y = "", x = 'Tissues in GTEx') +
  scale_fill_manual(values = count$ColorPlot) +
  theme(legend.position = "none")
p2


library(scales)
sample.size <- read.table('Data/00_sample_size_with_phenotype.txt', header = T)
count <- merge(count, sample.size, by = 'tissue')
count <- count[order(count$rhyGene.num),]
count$tissue <- factor(count$tissue, levels = count$tissue)
head(count)
median(count$samplesize)

ggplot(data = count, aes(x = samplesize, y = rhyGene.num)) + geom_point()
ggplot(data = count, aes(x = samplesize, y = rhyQTL.num)) + geom_point()

# heatmap of sample size
p3 <- ggplot(data = count, aes(x = 1, y = tissue, fill = samplesize)) + geom_tile() + 
  scale_fill_gradientn(colours = c("yellow", "orange", "darkred")) + theme_classic()

p3 

library(grid)
pdf('./13_Figures/00_Figure_rhyGene.rhyQTL.number.pdf', width = 14, height = 5.2, useDingbats = T)
print(p1, vp=viewport(.2, 1, x = .2, y = .5))
print(p2, vp=viewport(.2, 1, x = .5, y = .5))
print(p3, vp=viewport(.32, 1, x = .8, y = .5))
dev.off()

####################################################################################################
#
#                                        Figure 2b. tissue sharing rhyGenes
#
####################################################################################################
#rm(list = ls()) 
library(tidyverse)
library(data.table)
library(grid)

#######################
#
# classify rGenes
#
#######################
setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')

# classify rGenes
gene.list<- fread('13_Results/00_tissue_rhyGene.name.list')
common <- as.data.frame(table(gene.list$name))
head(common)
names(common) <- c('name', 'Tissue.num')

annot <- read.table('../gencode.v26.GRCh38.genes.annot')
annot <- unique(annot[, c('V6', "V7")])
names(annot) <- c('rhyGene', 'name')

common <- merge(common, annot, by = 'name')

common <- common[order(-common$Tissue.num),]

common$rhyGene.type <- ifelse(common$Tissue.num >= 20, 'Ubiquitous', 'Intermediately specific')
common$rhyGene.type[common$Tissue.num == 1] <- 'Tissue specific'

# number in fig2b
share.num <- as.data.frame(table(common$rhyGene.type))
share.num$ratio <- share.num$Freq/sum(share.num$Freq)
share.num


outdir <- paste0(dir, cat, '13_Results/02_tissue_sharing_rhyGene/')
if (!file.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
write.table(common, paste0(outdir, 'tissue_sharing_rhyGenes.txt'), quote = F, sep = '\t', row.names = F)
write.table(unique(filter(common, Tissue.num >= 20)$name), paste0(outdir, 'Ubiquitous.rhyGene.txt'), quote = F, sep = '\t', row.names = F, col.names = F)


### tissue specific rhyGenes
outdir <- paste0(dir, cat, '13_Results/02_tissue_sharing_rhyGene/tissue_specific_rhyGene/')
if (!file.exists(outdir)) {dir.create(outdir, recursive = TRUE)}

gene.list <- fread('13_Results/00_tissue_rhyGene.name.list')
for(x in unique(gene.list$tissue)){
  
  write.table(intersect(gene.list$name[gene.list$tissue == x], common$name[common$Tissue.num == 1]),
              paste0(outdir, x,'.unique.rhyGene.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
  
}

#################
#
# rhyGenes category heatmap
#
##################
# input is rhyGene or not?
gene.list <- fread('13_Results/00_tissue_rhyGene.name.list')
names(gene.list) <- c('rhyGene', 'name','tissue')
gene.list$yorn <- 1
gene.list.convert <- spread(gene.list, tissue, yorn)
gene.list.convert[is.na(gene.list.convert)] <- 0

# get data for plot 
r <- gather(gene.list.convert, tissue, yorn, `Adipose-Subcutaneous`:Vagina)


# rank tissue according to # of rhyGene
num <- as.data.frame(table(gene.list$tissue))
names(num) <- c('tissue', 'rhyGene.num')
num <- num[order(num$rhyGene.num),]
num$tissue <- factor(num$tissue, levels = num$tissue)
tissue_order <- num$tissue 
r$tissue <- factor(r$tissue, levels = tissue_order)


# rank rhyGene 
## 1. rank rhyGene in >=2 tissues
gene_order1 <- common$rhyGene[order(-common$Tissue.num[common$Tissue.num >= 2])]


# 2. rank tissue specific rhyGenes according to the tissue rand
common <- fread('13_Results/02_tissue_sharing_rhyGene/tissue_sharing_rhyGenes.txt')
tissue_specific <- filter(common, Tissue.num == 1)

# get gene order in only 1 tissue
gene_order <- c()
tissue_specific_rhyGene <- c()
for(x in rev(tissue_order)){
  
  tmp <- tissue_specific$rhyGene[match(intersect(gene.list$rhyGene[gene.list$tissue == x], tissue_specific$rhyGene), tissue_specific$rhyGene)]
  gene_order <- c(gene_order, tmp)
  tissue_specific_rhyGene <- data.frame(tissue = x, rhyGene = tmp) %>% rbind(tissue_specific_rhyGene, .)
  
}

# output tissue specific genes
write.table(tissue_specific_rhyGene, '13_Results/02_tissue_sharing_rhyGene/tissue_specific_rhyGene.txt', row.names = F, quote = F, sep = '\t')

r$rhyGene <- factor(r$rhyGene, levels = c(gene_order1, gene_order))
gene.rank <- data.frame(rhyGene = levels(r$rhyGene), rank = 1:length(levels(r$rhyGene)))
r <- merge(r, gene.rank, by = 'rhyGene')
head(r)

theme <-  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 8),
        axis.title = element_text(size = 8),
        plot.title = element_text(size = 8),
        axis.line = element_line(linewidth = 0.2),
        axis.ticks = element_line(linewidth = 0.2)) + 
  theme(legend.position = "none")

p <- ggplot(data = r, aes(x = rank, y = tissue, fill = as.character(yorn))) + 
  geom_tile() +   
  theme + 
  scale_fill_manual(values = c("grey100", "#FA8072")) +
  #theme(axis.text.x = element_blank())+
  #axis.ticks.x = element_blank())+
  labs(x = '', y = "") +
  geom_vline(xintercept = nrow(filter(common, Tissue.num >= 20)))+
  geom_vline(xintercept = nrow(filter(common, Tissue.num >= 2)))
p

plot.outdir <- paste0(dir, cat, '13_Figures/02_tissue_sharing_rhyGene/')
if (!file.exists(plot.outdir)) {dir.create(plot.outdir, recursive = TRUE)}
pdf(paste0(plot.outdir, 'tissue_sharing_rhyGene.pdf'), width = 5, height = 5.2, useDingbats = F)
print(p, vp=viewport(.9, .9, x = .5, y = .5))
dev.off()

## get tissue_specific_rhyGene.num
tissue_specific_rgene <- fread('13_Results/02_tissue_sharing_rhyGene/tissue_specific_rhyGene.txt')
a <- as.data.frame(table(tissue_specific_rgene$tissue))
names(a) <- c('tissue', 'tissue_specific_rhyGene.num')
head(a)

####################################################################################################
#
#                                               Figure 2c-g. phase distribution
#
####################################################################################################
dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/"
cat <- "cis_QTL/"

file <- list.dirs(paste0(dir, cat, '00_rQTL_mapping/00_Genotype'))
tissue <- gsub(paste0(dir, cat, '00_rQTL_mapping/00_Genotype/'), '', file)[-1]


outdir <- '13_Results/03_rhyGene_phase/rhyGene_phase_12h_interval/'
dir.create(outdir, recursive = TRUE)
dir.create('13_Figures/03_rhyGene_phase/', recursive = TRUE)

all <- c()
pdf(paste0('./13_Figures/03_rhyGene_phase/all_catagory_distribution.pdf'))
for(i in 1:length(tissue)){
  
  tmp <- fread(paste0(dir, cat, '13_Results/00_all_rhyQTL/', tissue[i], '.txt'))
  tmp$phase <- apply(select(tmp, contains("phase"), max_amp_idx), 1, function(x) (x[x[3]]))
  mean.phase <- as.data.frame(tmp %>% group_by(gene) %>% summarise(mean.phase = median(phase)))
  
  annot <- fread(paste0(dir, cat, 'all_gene.annot'))
  mean.phase <- merge(mean.phase, annot, by.x = 'gene', by.y = 'V6')
  mean.phase$tissue <- tissue[i]
  all <- rbind(all, mean.phase[,c(1,2,8,12)])
  
  
  write.table(filter(mean.phase, mean.phase > 0 & mean.phase <= 12)[,8], paste0(outdir, tissue[i], '.AM.rhyGene.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
  write.table(filter(mean.phase, mean.phase > 12 & mean.phase <= 24)[,8], paste0(outdir, tissue[i], '.PM.rhyGene.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
  
  print(ggplot(data = mean.phase, aes(x = mean.phase)) + geom_density() + theme_minimal() + ggtitle(tissue))
}
dev.off()  
write.table(all, '13_Results/03_rhyGene_phase/all.tissue.rhyGene.phase.txt', quote = F, sep = '\t', row.names = F, col.names = F)  


############## plot all tissue in one figure ###############################
all <- fread('13_Results/03_rhyGene_phase/all.tissue.rhyGene.phase.txt')
names(all) <- c('gene', 'mean.phase', 'name', 'tissue')
color <- fread('Data/color_annotation_v3.txt')

all$tissue <- factor(all$tissue, levels=unique(all$tissue))
a <- as.data.frame(table(all$tissue))
ggplot(data = filter(all, tissue != 'Brain-Substantianigra'), aes(x = mean.phase, color = tissue)) + geom_density() + coord_polar() + theme_minimal() +
  scale_colour_manual(values = color$ColorPlot[match(levels(all$tissue), color$tissue_rQTL)]) + 
  theme(legend.position = "none") +
  scale_x_continuous(breaks = seq(0, 24, 4))+
  theme(axis.text= element_text(color = "black", size = 8),
        axis.title = element_text(size = 8))

###################################################### PLOT all class of tissues ####################################
levels(all$tissue)
all <- fread('13_Results/03_rhyGene_phase/all.tissue.rhyGene.phase.txt')
names(all) <- c('gene', 'mean.phase', 'name', 'tissue')
color <- fread('Data/color_annotation_v3.txt')
all <- filter(all, tissue != 'Brain-Substantianigra')

outdir <- paste0(dir, cat, '13_Figures/03_rhyGene_phase/')
for(x in unique(color$Class)[unique(color$Class) !=  "Cells"]){
  
  tmp <- filter(all, tissue %in% color$tissue_rQTL[color$Class == x])
  tmp$tissue <- factor(tmp$tissue, levels=unique(tmp$tissue))
  
  pdf(paste0(outdir, 'all_class.', gsub('/', '', x),'_phase_distribution.pdf'), width = 5, height = 2, useDingbats = F)
  print(
    ggplot(data = tmp, aes(x = mean.phase, color = tissue)) + geom_density() + coord_polar() + theme_minimal() +
      scale_colour_manual(values = color$ColorPlot[match(levels(tmp$tissue), color$tissue_rQTL)]) + 
      scale_x_continuous(limits = c(0, 24), breaks = c(0, 4, 8, 12, 16, 20, 24)) + 
      labs(x = 'Time of day (h)', y = 'Density')+
      theme(axis.line = element_blank(),
            axis.text = element_text(color = "black", size = 8),
            axis.title = element_text(size = 8),
            plot.title = element_text(size = 8)) +
      ggtitle(x)
  )
  dev.off()
}

outdir <- paste0(dir, cat, '13_Figures/03_rhyGene_phase/')
for(x in unique(color$Class)[unique(color$Class) !=  "Cells"]){
  
  tmp <- filter(all, tissue %in% color$tissue_rQTL[color$Class == x])
  tmp$tissue <- factor(tmp$tissue, levels=unique(tmp$tissue))
  
  pdf(paste0(outdir, 'all_class.', gsub('/', '', x),'_phase_distribution_dashed.pdf'), width = 5, height = 2, useDingbats = F)
  print(
    ggplot(data = tmp, aes(x = mean.phase, color = tissue)) + geom_density(linetype = 'dashed') + coord_polar() + theme_minimal() +
      scale_colour_manual(values = color$ColorPlot[match(levels(tmp$tissue), color$tissue_rQTL)]) + 
      scale_x_continuous(limits = c(0, 24), breaks = c(0, 4, 8, 12, 16, 20, 24)) + 
      labs(x = 'Time of day (h)', y = 'Density')+
      theme(axis.line = element_blank(),
            axis.text = element_text(color = "black", size = 8),
            axis.title = element_text(size = 8),
            plot.title = element_text(size = 8)) +
      ggtitle(x)
  )
  dev.off()
}

############
# annotation tissue name with color
############
color <- fread('Data/color_annotation_v3.txt')
for(x in unique(color$Class)[unique(color$Class) !=  "Cells"]){
  
  tmp <- filter(color, Class == x)
  #tmp$FullName <- factor(tmp$FullName, levels = tmp$FullName)
  tmp <- tmp[order(rev(tmp$FullName)),]
  tmp$x <- 1
  tmp$y <- 1:nrow(tmp)
  tmp$FullName <- factor(tmp$FullName, levels = tmp$FullName)
  tmp
  
  pdf(paste0('./13_Figures/03_rhyGene_phase/all.class.', gsub('/', '', x),'_annot_font.pdf'), width = 5, height = 2, useDingbats = F)
  print(
    
    ggplot(data = tmp, aes(x = x, y = y)) + geom_point() +
      scale_x_continuous(limits = c(0, 24), breaks = c(0, 6, 12, 18, 24)) +
      theme(axis.line = element_blank(),
            axis.text = element_text(color = "black", size = 8),
            #axis.text = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            axis.title = element_text(size = 8),
            plot.title = element_text(size = 8)) +
      annotate("text", x = (tmp$x + 0.1), y = tmp$y, label = tmp$FullName, color = tmp$ColorPlot[match(levels(tmp$FullName), tmp$FullName)], size = 3) +
      theme_void() +
      ggtitle(x)
    
  )
  dev.off()
}



################################## plot phase distribution in liver tissue ###################################################
pdf(paste0('./13_Figures/03_rhyGene_phase/Liver.pdf'), width = 2, height = 2, useDingbats = F)
print(
  ggplot(data = filter(all, tissue == 'Liver'), aes(x = mean.phase)) + geom_density(color = color$ColorPlot[match('Liver', color$tissue_rQTL)]) + 
    coord_polar() + theme_minimal() +
    scale_x_continuous(limits = c(0, 24), breaks = c(0, 4, 8, 12, 16, 20, 24)) + 
    labs(x = 'Time of day (h)', y = 'Density')+
    theme(axis.line = element_blank(),
          axis.text = element_text(color = "black", size = 8),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 8)) 
)
print(
  ggplot(data = filter(all, tissue == 'Liver'), aes(x = mean.phase)) + geom_density(color = color$ColorPlot[match('Liver', color$tissue_rQTL)], linetype = 'dashed') + 
    coord_polar() + theme_minimal() +
    scale_x_continuous(limits = c(0, 24), breaks = c(0, 4, 8, 12, 16, 20, 24)) + 
    labs(x = 'Time of day (h)', y = 'Density')+
    theme(axis.line = element_blank(),
          axis.text = element_text(color = "black", size = 8),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 8)) 
)
dev.off()



#########################################  skin tissues ########################################

# skin AM rhyGene compare
m.1 <- fread('13_Results/03_rhyGene_phase/rhyGene_phase_12h_interval/Skin-SunExposed_Lowerleg.AM.rhyGene.txt', header = F)
m.2 <- fread('13_Results/03_rhyGene_phase/rhyGene_phase_12h_interval/Skin-NotSunExposed_Suprapubic.AM.rhyGene.txt', header = F)

m.common.am <- intersect(m.1$V1, m.2$V1)
sun.unique.am <- setdiff(m.1$V1, m.common.am)
no.sun.unique.am <- setdiff(m.2$V1, m.common.am)

dir.create('13_Results/03_rhyGene_phase/skin/', recursive = TRUE)
write.table(m.common.am, '13_Results/03_rhyGene_phase/skin/Skin-AM-common-rhyGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(sun.unique.am, '13_Results/03_rhyGene_phase/skin/Skin-AM-SunExposed-unique-rhyGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(no.sun.unique.am, '13_Results/03_rhyGene_phase/skin/Skin-AM-Nosun-unique-common-rhyGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)

length(m.common.am)
length(sun.unique.am)
length(no.sun.unique.am)

# skin PM rhyGene compare
m.1 <- fread('13_Results/03_rhyGene_phase/rhyGene_phase_12h_interval/Skin-SunExposed_Lowerleg.PM.rhyGene.txt', header = F)
m.2 <- fread('13_Results/03_rhyGene_phase/rhyGene_phase_12h_interval/Skin-NotSunExposed_Suprapubic.PM.rhyGene.txt', header = F)

m.common.pm <- intersect(m.1$V1, m.2$V1)
sun.unique.pm <- setdiff(m.1$V1, m.common.pm)
no.sun.unique.pm <- setdiff(m.2$V1, m.common.pm)

length(m.common.pm)
length(sun.unique.pm)
length(no.sun.unique.pm)

write.table(m.common.pm, '13_Results/03_rhyGene_phase/skin/Skin-PM-common-rhyGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(sun.unique.pm, '13_Results/03_rhyGene_phase/skin/Skin-PM-SunExposed-unique-rhyGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(no.sun.unique.pm, '13_Results/03_rhyGene_phase/skin/Skin-PM-Nosun-unique-common-rhyGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)
