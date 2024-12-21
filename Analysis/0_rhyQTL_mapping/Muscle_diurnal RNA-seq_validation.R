rm(list=ls())
library(data.table)
library(tidyverse)


#############################################
#
# 0. process vcf file AA -> 0 AC -> 1 CC -> 2
#
#############################################
vcf <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/04_validation_muscle/Data/VCF/snp_filtered.vcf.gz') %>% as.data.frame()

# only keep genotype information
sub <- select(vcf, contains('ZT12'))
result <- as.data.frame(apply(sub, c(1, 2), function(row) {
  str_split_fixed(row, ':', 5)[1,1]
}))

# get genotype for each individual
split_and_sum <- function(element) {
  element <- gsub('\\.', '9', element)
  split_elements <- strsplit(element, "[|/]")  
  sum_of_elements <- sapply(split_elements, function(x) sum(as.numeric(x)))
  return(sum_of_elements)
}
result_sum <- as.data.frame(t(apply(result, 1, split_and_sum)))

vcf.new <- cbind(vcf[,c(1,2,4,5)], result_sum)

names(vcf.new)[grep("CHROM", names(vcf.new))] <- 'CHR'

vcf.new$ID <- paste(vcf.new$CHR, vcf.new$POS, vcf.new$REF, vcf.new$ALT, sep = '_')
vcf.new <- select(vcf.new, ID, CHR, POS, REF, ALT, contains('ZT'))

write.table(vcf.new, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/04_validation_muscle/Data/VCF/snp_processed.vcf', quote = F, sep = '\t', row.names = F)

######################################################
#
# 1. generate genotype information for each SNP 
#
######################################################
vcf.new <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/04_validation_muscle/Data/VCF/snp_processed.vcf') %>% as.data.frame()
qtl <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/00_all_rhyQTL/Muscle-Skeletal.txt')

length(intersect(gsub("_b38", "", qtl$ID), vcf.new$ID))

# only keep rQTLs in muscle
vcf <- filter(vcf.new, ID %in% intersect(gsub("_b38", "", qtl$ID), vcf.new$ID)) 

geno <- lapply(vcf$ID, function(x){
  
  tmp.vcf <- as.data.frame(t(filter(vcf, ID == x)[,c(6:ncol(vcf))]))
  names(tmp.vcf) <- 'geno'
  return(filter(tmp.vcf, geno != 18) %>% rownames_to_column(., var = 'sample'))
  
})
names(geno) <- vcf$ID
saveRDS(geno, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/04_validation_muscle/Data/geno.rds')

##############################################################
#
# 2. count # of individuals in each genotype for each SNP id 
#
##############################################################

geno <- readRDS('04_validation_muscle/Data/geno.rds')
geno.count <- lapply(geno, function(x){
  
  num <- as.numeric(table(factor(x$geno, levels = c(0, 1, 2))))
  count <- data.frame(geno_0 = num[1],
                      geno_1 = num[2],
                      geno_2 = num[3])
  return(count)
  
})
geno.count <- as.data.frame(do.call(rbind, geno.count)) %>% rownames_to_column(., var = 'ID')


###################################################
#
# 3. filter SNP id top 2 genotype with sample number >= 2
# 
###################################################
id.filter <- lapply(geno.count$ID, function(x){
  
  if(sort(as.numeric(filter(geno.count, ID == x)[,2:4]))[2] >= 2){
    
    tmp <- filter(geno.count, ID == x) %>%
      mutate(geno.no.1 = sort(as.numeric(filter(geno.count, ID == x)[,2:4]))[3]) %>%
      mutate(geno.no.2 = sort(as.numeric(filter(geno.count, ID == x)[,2:4]))[2]) %>%
      mutate(geno.no.3 = sort(as.numeric(filter(geno.count, ID == x)[,2:4]))[1])
    
    
    return(tmp)
    
  }else{
    
    return(NULL)
  }
  
})

id.filter <- as.data.frame(do.call(rbind, id.filter))
write.table(id.filter, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/04_validation_muscle/Data/VCF/snp_id.filter_gt2.txt', quote = F, sep = '\t', row.names = F)


#######################################################
#
# 4. get top 2 genotypes
#
#######################################################
id.filter <- fread('04_validation_muscle/Data/VCF/snp_id.filter_gt2.txt')
geno <- readRDS('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/04_validation_muscle/Data/geno.rds')

top_two <- apply(id.filter[,2:4], 1, function(row) {
  sorted_row <- names(sort(row, decreasing = TRUE))[1:2]
  return(sorted_row)
})
top_two <- as.data.frame(t(top_two))
top_two$ID <- id.filter$ID
names(top_two) <- c('top1_geno', 'top2_geno', 'ID')
top_two <- merge(id.filter, top_two, by = 'ID')
write.table(top_two, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/04_validation_muscle/Data/VCF/snp_id.filter_gt2_geno_info.txt', quote = F, sep = '\t', row.names = F)


geno <- readRDS('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/04_validation_muscle/Data/geno.rds')
top_two <- read.table('04_validation_muscle/Data/VCF/snp_id.filter_gt2_geno_info.txt', header = T)
###############################
#
# get top 2 genotype for rhyGenes-rhyQTL pairs
#
################################
expr <- fread('04_validation_muscle/Data/GSE108539_norm_reads_count.txt.gz') %>% as.data.frame() %>% select(Gene_name, contains('Exon')) 
rownames(expr) <- expr$Gene_name

qtl <- fread('13_Results/00_all_rhyQTL/Muscle-Skeletal.txt')

cir.compare <- lapply(top_two$ID, function(tmp.id){
  
  print(tmp.id)
  top1 <- gsub('geno_', '', top_two$top1_geno[top_two$ID == tmp.id])
  top2 <- gsub('geno_', '', top_two$top2_geno[top_two$ID == tmp.id])
  
  tmp.sample <- filter(geno[[tmp.id]], geno == top1 | geno == top2)
  tmp.sample$sample <- gsub('_ZT12', '', tmp.sample$sample) %>% gsub("^", "Exon_", .) %>% gsub("S[0-9]+", "", .)
  
  data <- select(expr, contains(tmp.sample$sample))
  
  group <- data.frame(sample = paste0(str_split_fixed(names(data), '_', 4)[ ,1], "_",
                                      str_split_fixed(names(data), '_', 4)[ ,2], "_"), 
                      id = names(data))
  group <- merge(group, tmp.sample, by = 'sample')
  group$time <- as.numeric(str_split_fixed(group$id, '_', 4)[ ,4])
  
  gene <- str_split_fixed(qtl$gene[qtl$ID == paste0(tmp.id, "_b38")], '\\.', 2)[ ,1]
  
  group.list <- lapply(gene, function(gene.tmp){
    
    if(sum(grepl(gene.tmp, rownames(expr))) == 0){
      
      return(NULL)
      
    }else{
      
      data <- as.data.frame(t(select(expr, contains(group$sample))[gene.tmp,])) %>% rownames_to_column(., var = 'id')
      names(data) <- c('id', 'expr')
      
      group <- merge(group, data, by = 'id')
      group$geno <- as.character(group$geno)
      
      return(group)
    }
    
  })
  
  names(group.list) <- gene
  group.list <- as.data.frame(do.call(rbind, group.list)) %>% rownames_to_column(., var = 'gene')
  group.list$gene <- gsub('\\.[0-9]+', '', group.list$gene)
  return(group.list)
  
})

names(cir.compare) <- top_two$ID
saveRDS(cir.compare, '04_validation_muscle/muscle.rhyGene.top2.genotype.compare.circadian.RDS')



##############################
#
#  fit circadian curve
#
############################

###################### Function ###################################
# Harmonic regression function
harm_reg <- function(x, t, period){
  n=length(x)
  fit0=lm(x~1)
  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)
  fit1=lm(x~c+s)
  mu=coef(fit1)[1]
  a=coef(fit1)[2]
  b=coef(fit1)[3]
  p.val=lrtest(fit1, fit0)$Pr[2]
  amp=2*sqrt(a^2+b^2)
  phase=atan2(b,a)%%(2*pi)
  phase=period*phase/(2*pi)
  return(c(pval=p.val,phase=phase,amp=amp,mu=mu,a=a,b=b, period=period, R2=summary(fit1)$r.squared))
}

suppressMessages(library("lmtest"))

#cir.compare <- readRDS('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/04_validation_muscle/Data/cir_data.rds')
cir.compare <- compact(cir.compare)

period <- 24
para <- lapply(cir.compare, function(tmp.2){
  
  para.list <- c()  
  for(gene.tmp in unique(tmp.2$gene)){
    
    tmp <- filter(tmp.2, gene == gene.tmp)
    regression <- as.numeric(unique(tmp$geno))
    
    for(genotype in unique(tmp$geno)){
      
      tmp.1 <- filter(tmp, geno == genotype)
      x <- tmp.1$expr
      t <- tmp.1$time
      
      regression <- c(regression, as.numeric(harm_reg(x, t, period)))
    }
    
    regression <- as.data.frame(t(regression))
    names(regression) <- c('geno_1', 'geno_2',
                           paste0(c("pval","phase.s","amp.c","mu.(Intercept)","a.c","b.s","period","R2"), "_0"),
                           paste0(c("pval","phase.s","amp.c","mu.(Intercept)","a.c","b.s","period","R2"), "_1"))
    
    regression$gene <- gene.tmp
    
    para.list <- rbind(para.list, regression)
  }
  
  
  return(para.list)
})

para <- as.data.frame(do.call(rbind, para)) %>% rownames_to_column(., var = 'id')
para$id <- gsub("\\.[0-9]+", "", para$id)
length(unique(para$gene))
para$cir.num <- apply(select(para, contains('pval')), 1, function(row) sum(row < 0.05))
length(unique(filter(para, cir.num != 0)$gene))

write.table(para, '04_validation_muscle/muscle.rhyGene.regression_in_top2_genotype.txt', quote = F, sep = '\t', row.names = F)


#####################################################
#
# 7. get percentage of validated rhyGene/rhyQTLs
#
######################################################
para <- fread('04_validation_muscle/muscle.rhyGene.regression_in_top2_genotype.txt')
top_two <- read.table('04_validation_muscle/Data/VCF/snp_id.filter_gt2_geno_info.txt', header = T)
names(top_two) <- c('ID', 'geno_0_num', 'geno_1_num', 'geno_2_num', names(top_two)[5:9])
para <- merge(top_two, para, by.x = 'ID', by.y = 'id')
para$phase_dif <- abs(para$phase.s_0 - para$phase.s_1)
para$amp_dif <- abs(para$amp.c_0 - para$amp.c_1)

hist(para$phase_dif)
para$type <- 'no_circadian'
para$type[para$cir.num == 1] <- '1_circadian'
para$type[para$cir.num == 2] <-  '2_circadian_no_diff'
para$type[para$cir.num == 2 & para$amp_dif >= 0.58] <- '2_circadian_diff'
para$type[para$cir.num == 2 & para$phase_dif >= 3] <- '2_circadian_diff'
para$cir.compare <- para$type
para$cir.compare[para$cir.compare == '2_circadian_diff'|para$cir.compare == '1_circadian'] <- 'diff_circadian'

table(para$type)
table(para$cir.compare)

# individual rhythm fit
rhy.convert <- fread('04_validation_muscle/Data/rhythm_fit_intron_exon.txt')
filter <- filter(rhy.convert, p.val_sig == 2 & amp_gt1.5 >= 2 & sample.id != '_7_') 
rhy.fit.num <- as.data.frame(table(filter$gene))
head(rhy.fit.num )
length(rhy.fit.num$Var1[rhy.fit.num$Freq > 1])

####################### use top 2 genotype gt 2
para.2 <- filter(para, gene %in% unique(rhy.fit.num$Var1[rhy.fit.num$Freq > 1]))
para.2.validate <- filter(para.2, cir.compare == 'diff_circadian')


table(para.2$cir.compare)

# # of ID and gene detected 
# original 1387, 292
# new  2366 229
length(unique(para.2$ID))
length(unique(para.2$gene))

# ratio of no circadian in both genotype
length(unique(para.2$ID[para.2$cir.compare == 'no_circadian']))/length(unique(para.2$ID))

# # of ID and gene validated
# original 754 , 225
# new 1341 193
# 
length(unique(para.2.validate$ID))
length(unique(para.2.validate$gene))

# validated gene ratio
length(unique(para.2.validate$gene))/length(unique(para.2$gene))
# validated rhyQTL ratio
length(unique(para.2.validate$ID))/length(unique(para.2$ID))

