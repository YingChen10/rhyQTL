# annot <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/all_gene.annot')
# annot <- annot[,c(1:5,7)]
# names(annot) <- c('chr', 'start', 'end', 'sense', 'gene', 'name')
# annot$tss <- annot$start
# annot$tss[annot$sense == '-'] <- annot$end[annot$sense == '-']
# 
# bed <- data.frame(c = annot$chr,
#                   s = annot$tss - 5000,
#                   e = annot$tss + 5000,
#                   n = annot$name)
# 
# bed <- filter(bed, c != 'chrX' & c != 'chrY' & c != 'chrM')
# write.table(bed, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/05_trans_rGene/gene_5k_flanking.bed', quote = F, sep = '\t', row.names = F, col.names = F)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# library(BiocManager)
# BiocManager::install("RCircos")

setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/08_trans_rhyQTL')

library(data.table)
library(tidyverse)
tissue <- list.dirs('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/08_trans_rhyQTL') %>% 
  gsub("/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/08_trans_rhyQTL/", '',.)

tissue <- tissue[2:length(tissue)]

# combine all trans-rhyGenes
all <- c()
for(x in tissue){
  
  trans <- fread(paste0(x, "/NR1D1.trans_gene.txt"))
  trans$tissue <- x
  all <- rbind(all, trans)
  
}

all$max_amp <- apply(select(all, contains("amp.c")), 1, function(x) max(x))
all$max_amp_idx <- apply(select(all, contains("amp.c_")), 1, function(x) which.max(x))
all$pval <- apply(select(all, contains("pval_"), max_amp_idx), 1, function(x) (x[x[3]]))
all$qval <- apply(select(all, contains("qval_"), max_amp_idx), 1, function(x) (x[x[3]]))
all$HANOVA.norm.qval <- p.adjust(all$HANOVA.norm, "BH")

all <- filter(all, pval <= 5e-4 & HANOVA.norm.qval < 0.01)
length(unique(all$name))
all <- filter(all, name != 'NR1D1')

#write.table(a, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/08_trans_rhyQTL/Brain.NR1D1.trans_gene.filter.txt', quote = F, sep = '\t', row.names = F)

length(unique(all$name))

all$idx <- paste0(all$tissue, '_', all$ID, "_", all$gene)
all.1 <- all
# remove cir-rhyGenes
for(x  in tissue){
  
  rqtl <- fread(paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/00_all_rhyQTL/', x, '.txt'))
  index <- paste0(x, '_', rqtl$ID, '_', rqtl$gene)
  all <- filter(all, !(idx %in% index))
  
}

length(unique(all$gene))

write.table(all, 'NR1D1_trans_rhyGene_brain.txt', quote = F, sep = '\t', row.names = F)

######### prepare for plot
all <- fread('NR1D1_trans_rhyGene_brain.txt')
head(all)

tables2 <- unique(select(all, gene, name))
write.table(tables2, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/08_trans_rhyQTL/TableS2.NR1D1_trans_rhyGene_brain.txt', quote = F, sep = '\t', row.names = F)

# generate circle plot line data 
line <- select(all, ID, gene)
pos <- read.table('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/all_gene.annot')
pos <- pos[,c(5, 1, 2, 3)]
names(pos) <- c('gene', 'Chromosome.1', 'chromStart.1', 'chromEnd.1')

line$Chromosome <- str_split_fixed(line$ID, '_', 5)[,1]
line$chromStart <- as.numeric(str_split_fixed(line$ID, '_', 5)[,2])
line$chromEnd <- as.numeric(str_split_fixed(line$ID, '_', 5)[,2])
line <- merge(line[,c(2:5)], pos, by = 'gene')

line <- line[,2:7]
line$chromEnd <- line$chromStart + 500
line$chromEnd.1 <- line$chromStart.1 + 500
#line <- column_to_rownames(line, var = 'gene')



# gene annotation
gene <- read.table('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/human_coding_gene_TSS_flanking1M.txt', header = T)
gene <- filter(gene, gene_id %in% c(all$gene, 'ENSG00000126368.5'))
gene <- select(gene, chr, start, end, gene_name)
names(gene) <- c('Chromosome', 'chromStart', 'chromEnd', 'Gene')
gene$Gene <- as.factor(gene$Gene)
gene$Chromosome <- as.factor(gene$Chromosome)
gene <- filter(gene, Gene %in% c('NR1D1', 'PER1', 'PER2', 'PER3', 'CRY1', 'ARNTL'))


library(RCircos)
# https://cloud.tencent.com/developer/article/1667791
data("UCSC.HG38.Human.CytoBandIdeogram")
cyto.info <- UCSC.HG38.Human.CytoBandIdeogram

# excluded chromosome
chr.exclude <- NULL
# # of circle
tracks.inside <- 10
#  
tracks.outside <- 0

RCircos.Set.Core.Components(cyto.info, chr.exclude, tracks.inside, tracks.outside)  
RCircos.List.Plot.Parameters()

RCircos.Set.Plot.Area()    
RCircos.Chromosome.Ideogram.Plot()

# add gene annotation line
side <- "in"
track.num <- 1
RCircos.Gene.Connector.Plot(gene, track.num, side)

# add gene name in the circle
name.col <- 4
track.num <- 2

# plot
RCircos.Gene.Name.Plot(gene, name.col, track.num, side)

# add line
track.num <- 5
RCircos.Link.Plot(line, track.num, TRUE)

