library(ggrepel)
library(pheatmap)
library(grid)
library(gridExtra)

theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) 


setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')
or <- fread('Result/02_Functional_annotation_enrichment/motif_enrichment_OR/homer.or.tissue.txt')

###################
#
# plot heatmap
#
##################
or$TF <- or$name
or$TF <- gsub("[^A-Za-z0-9-].*", "", or$TF)

or$rql.num <- or$rqtl.motif + or$rqtl.non.motif
or.filter <- filter(or, rql.num > 2000)
num <- filter(or, TF == 'RUNX1')


# E-box: get max CACGTG odds ratio 
filter <- or[grep('CACGTG', or$motif),]

ebox <- as.data.frame(
  filter %>%
    group_by(tissue) %>%
    slice(which.max(odds.ratio)))
ebox$motif <- 'CACGTG'


# get reverb odds ratio
reverb <- filter(or, TF == 'Reverb')


# get HSF
filter.1 <- or[grep('TTCTAGAA', or$motif),]
hsf <- as.data.frame(
  filter.1 %>%
    group_by(tissue) %>%
    slice(which.max(odds.ratio)))


combine <- rbind(ebox, reverb, hsf)
color <- fread('Data/color_annotation_v3.txt')
combine <- merge(combine, color, by.x = 'tissue', by.y = 'tissue_rQTL')
write.table(combine, 'Result/02_Functional_annotation_enrichment/motif_enrichment_OR/homer.or.tissue.heatmap_input.txt', quote = F, sep = '\t', row.names = F)

combine <- fread('Result/02_Functional_annotation_enrichment/motif_enrichment_OR/homer.or.tissue.heatmap_input.txt')
head(combine)
combine <- combine %>% filter(rql.num > 2000) %>% select(motif, ShortName, odds.ratio)
combine$ShortName <- factor(combine$ShortName, levels = unique(combine$ShortName))

filter_convert <- combine %>% select(., motif, ShortName, odds.ratio) %>% spread(., motif, odds.ratio) 
rownames(filter_convert) <- filter_convert$ShortName
filter_convert <- filter_convert[,-1]


p1 <- pheatmap(filter_convert,  breaks = seq(0, 2.0, length.out = 101), cluster_rows = FALSE, cluster_cols = FALSE)
p1

pdf("Figure/02_Functional_annotation_enrichment/core_clock_heatmap.pdf", height  = 7, width = 2.5)
p1 <- pheatmap(filter_convert, breaks = seq(0, 2.0, length.out = 101), cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()


mean <- as.data.frame(combine %>% group_by(motif) %>% summarise(mean = mean(odds.ratio), sd = sd(odds.ratio)))
p2 <- ggplot(data = mean, aes(x = motif, y = mean)) + geom_col(fill = '#4682B4', width = 0.2) + geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0) + theme +
  geom_hline(yintercept = 1, color = 'grey30', linetype = 'dashed') + labs(y = 'Mean of odds ratios', x = '')
p2

library(grid)
pdf('./Figure/02_Functional_annotation_enrichment/homer_enriched_motif_or.pdf', width = 7, height = 4.5, useDingbats = F)
print(p2, vp=viewport(.2, .4, x =.6, y = .5))
dev.off()


###################
#
# plot tissue Adipose-Visceral_Omentum as an example
#
##################
or <- fread('Result/02_Functional_annotation_enrichment/motif_enrichment_OR/homer.or.tissue.txt')
or$TF <- or$name
or$TF <- gsub("[^A-Za-z0-9-].*", "", or$TF)
x <- 'Adipose-Visceral_Omentum'
a <- filter(or, tissue == x )

# no ebox 
no.ebox <- a[!grep('CACGTG', a$motif),]
# ebox
filter <- a[grep('CACGTG', a$motif),]
ebox <- filter[order(-filter$odds.ratio),][1]
ebox$TF <- 'BMAL1'

a <- rbind(ebox, no.ebox)
a <- a[order(a$p.value),]
a$rank.pval <- c(1:nrow(a))
a <- a[order(-a$odds.ratio),]
a$rank <- c(1:nrow(a))


library(ggrepel)
p2 <- ggplot() + 
  geom_point(data = a, aes(x = rank, y = odds.ratio), size = 0.1) +
  geom_point(data =  filter(a, TF %in% c('Reverb', 'BMAL1', 'Sp1')), aes(x = rank, y = odds.ratio), size = 0.1, color = 'red')+
  # geom_text_repel(data = filter(a, TF %in% c('Reverb', "RORa", 'NPAS2', 'CLOCK', 'BMAL1')), aes(x = rank, y = odds.ratio, label = TF), size = 3) + 
  theme +
  ggtitle(x) + labs(x = 'Motif rank', y = 'Enrichment of rQTLs') +
  geom_hline(yintercept = 1) 
p2


pdf('./Figure/02_Functional_annotation_enrichment/homer_enriched_motif_or_Adipose.pdf', width = 4, height = 4, useDingbats = F)
print(p3, vp=viewport(.45, .8, x = .5, y = .5))
dev.off()

