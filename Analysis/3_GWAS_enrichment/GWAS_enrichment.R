
library(tidyverse)
library(data.table)
library(parallel)
setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')
dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/"
cat <- "cis_QTL/"


# GWAS_Catalog
GWAS <- fread('03_enrich_GWAS/gwas_catalog_v1.0.2-associations_e111_r2024-01-19.tsv', header = T, sep = "\t", na.strings = c("", "NA", "N/A"), quote="")
# 569163
length(unique(GWAS$`DISEASE/TRAIT`))
# 29509

#############################################
#
#            filter GWAS SNP
#
#################################################
# 1. Only keep European population
idx <- grep("European", as.character(GWAS$`INITIAL SAMPLE SIZE`), ignore.case = TRUE)
GWAS <- GWAS[idx, ]
# 393414
length(unique(GWAS$`DISEASE/TRAIT`))
# 18129

# 2. remove insignificant SNPs (P value > 5*10-8)
GWAS$`P-VALUE` <- as.numeric(GWAS$`P-VALUE`)
GWAS <- filter(GWAS, `P-VALUE` < 5e-8)
# 318783
length(unique(GWAS.filter$`DISEASE/TRAIT`))
# 16348

# 3. remove SNPs in the human leukocyte antigen locus (hg38: chr6:29,723,339-33,087,199)

GWAS <- transmute(GWAS,
                  SNP=SNPS,
                  CHR_ID = CHR_ID,
                  CHR_POS = CHR_POS,
                  Trait=`DISEASE/TRAIT`,
                  Maped_trait=MAPPED_TRAIT,
                  Maped_trait_URI=MAPPED_TRAIT_URI,
                  RISK_ALLELE=`STRONGEST SNP-RISK ALLELE`,
                  OR_or_BETA=`OR or BETA`,
                  P_value=`P-VALUE`,
                  Reportd_gene=`REPORTED GENE(S)`,
                  PUBMEDID,
                  Study_accession = `STUDY ACCESSION`,
                  LINK)

# GWAS$RISK_ALLELE <- gsub("^\\w+\\-","",GWAS$RISK_ALLELE, perl=T)
# GWAS$RISK_ALLELE <- gsub("\\-","",GWAS$RISK_ALLELE, perl=T)
# GWAS$RISK_ALLELE[!grepl("A|G|C|T", GWAS$RISK_ALLELE)] <- NA
# GWAS$P_value <- as.numeric(GWAS$P_value)
# GWAS$OR_or_BETA <- as.numeric(GWAS$OR_or_BETA)


# # format some tags with multiple SNP in a row
# multiple <- GWAS[grep(";", GWAS$SNP),]
# split.list <- lapply(1:nrow(multiple), function (i){
#   
#   snp <- as.character(unlist(strsplit(as.character(gsub(" ", '', multiple[i, 1])), ";")))
#   chr <- as.character(unlist(strsplit(as.character(gsub(" ", '', multiple[i, 2])), ";")))
#   pos <- as.character(unlist(strsplit(as.character(gsub(" ", '', multiple[i, 3])), ";")))
#   risk <-  as.character(unlist(strsplit(as.character(gsub(" ", '', multiple[i, 7])), ";")))
#   gene <-  as.character(unlist(strsplit(as.character(gsub(" ", '', multiple[i, 10])), ",")))
#   
#   tmp <- c()
#   for(j in 1:length(snp)){
#     
#     tmp <- data.frame(SNP = snp[j],
#                       CHR_ID = chr[j],
#                       CHR_POS = pos[j],
#                       Trait = as.character(multiple[i, 4]),
#                       Maped_trait = as.character(multiple[i, 5]),
#                       Maped_trait_URI = as.character(multiple[i, 6]),
#                       RISK_ALLELE = risk[j],
#                       OR_or_BETA = as.numeric(multiple[i, 8]),
#                       P_value = as.numeric(multiple[i, 9]),
#                       Reportd_gene = gene[j],
#                       PUBMEDID = multiple[i, 11],
#                       Study_accession = multiple[i, 12],
#                       LINK = multiple[i, 13]) %>% rbind(tmp, .)
#   }
#   
#   return(tmp)
#   
# })
# 
# split <- as.data.frame(do.call(rbind, split.list))

# GWAS <- GWAS[!(grep(";", GWAS$SNP)),]
# GWAS.filter <- rbind(GWAS.1, split)

GWAS.1 <- GWAS[!(grep(";", GWAS$SNP)),]
GWAS.2 <- GWAS.1[!(grep(",", GWAS.1$SNP)),]
GWAS.3 <- GWAS.2[!(grep("x", GWAS.2$SNP)),]

GWAS.3$RISK_ALLELE <- str_split_fixed(GWAS.3$RISK_ALLELE, '-', 2)[ ,2]
GWAS.3 <- GWAS.3[!is.na(GWAS.3$CHR_ID), ]
GWAS.3 <- GWAS.3[!is.na(GWAS.3$CHR_POS), ]
GWAS.4 <- GWAS.3[grepl("^[ATCG]+$", GWAS.3$RISK_ALLELE), ]
GWAS.filter <- GWAS.4

# remove SNPs in the human leukocyte antigen locus (hg38: chr6:29,723,339-33,087,199)
GWAS.filter$CHR_POS <- as.numeric(GWAS.filter$CHR_POS)
GWAS.filter <- filter(GWAS.filter, !(CHR_ID == 6 & CHR_POS > 29723339 & CHR_POS < 33087199))
GWAS.filter$Trait_mod <- gsub(' ', '-', GWAS.filter$Trait) %>% gsub("[^A-Za-z0-9._-]", "", .)
write.table(GWAS.filter, './03_enrich_GWAS/GWAS_enrichment/GWAS_r2024-01-19.filter.txt', quote = F, sep = '\t', row.names = F)

length(unique(GWAS.filter$Study_accession))
#[1] 14684
length(unique(GWAS.filter$PUBMEDID))
#[1] 2295




#######################################################
#
#   count the number of traits overlapped with rhyQTLs
#
#######################################################

file <- list.files('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/00_rQTL_mapping/04_Results/cis_rhyQTL_tissue/', pattern = '.rQTL')

# filtered GWAS
gwas.filter <- fread('./03_enrich_GWAS/GWAS_enrichment/GWAS_r2024-01-19.filter.txt') %>% as.data.frame(.)

overlap <- lapply(file, function(x){
  
  rqtl <- fread(paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/00_rQTL_mapping/04_Results/cis_rhyQTL_tissue/', x))
  
  tmp <- filter(gwas.filter, SNP %in% rqtl$rsid)
  
  overlaped.trait.num <- as.data.frame(tmp %>% group_by(Study_accession) %>% summarise(count = length(SNP)))
  
  result <-  data.frame(
    
    tissue = gsub('.rQTL', '', x),
    Study_accession = unique(tmp$Study_accession)
    
  )
  
  result <- merge(result, overlaped.trait.num, by = 'Study_accession')
  
  return(result)
  
})
overlap <- as.data.frame(do.call(rbind, overlap))

overlap.filter <- filter(overlap, count >= 3)


# filtered GWAS with unique trait number
gwas.filter <- fread('./03_enrich_GWAS/GWAS_enrichment/GWAS_r2024-01-19.filter.txt') %>% as.data.frame(.)

filter <- filter(gwas.filter, Study_accession %in% overlap.filter$Study_accession) %>% select(., Trait_mod, Study_accession, Maped_trait_URI) %>% unique(.)

filter <- filter %>%
  distinct(Trait_mod, .keep_all = TRUE)

overlap.filter <- merge(overlap.filter, filter, by = 'Study_accession')
write.table(overlap.filter, 'Result/03_enrich_GWAS/rhyQTL_associated_traits_id.txt', quote = F, row.names = F, sep  = '\t')

length(unique(overlap.filter$Maped_trait_URI))
num <- overlap.filter %>% group_by(tissue) %>% summarise(trait.num = n_distinct(Trait_mod)) %>% as.data.frame()
write.table(num, 'Result/03_enrich_GWAS/rhyQTL_associated_traits.txt', quote = F, row.names = F, sep  = '\t')


## Get rsid 
gwas.filter <- fread('./03_enrich_GWAS/GWAS_enrichment/GWAS_r2024-01-19.filter.txt') %>% as.data.frame(.)
overlap.filter <- fread('Result/03_enrich_GWAS/rhyQTL_associated_traits_id.txt')
overlap.filter$Trait <- 1
overlap.filter$SNP <- 1

for(i in 1:nrow(overlap.filter)){
  
  print(i)
  tissue <- overlap.filter$tissue[i]
  rqtl <- fread(paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/00_rQTL_mapping/04_Results/cis_rhyQTL_tissue/', tissue, '.rQTL'))
  tmp <- filter(gwas.filter, SNP %in% rqtl$rsid)
  Study_accession.tmp <- overlap.filter$Study_accession[i]
  overlap.filter$Trait[i] <- filter(tmp, Study_accession == Study_accession.tmp)$Trait_mod[1]
  overlap.filter$SNP[i] <- paste(filter(tmp, Study_accession == Study_accession.tmp)$SNP, collapse = ", ")
  
}

write.table(overlap.filter, 'Result/03_enrich_GWAS/TableS6_ori_20240619.txt', sep = '\t', quote = F, row.names = F)



num <- fread('Result/03_enrich_GWAS/rhyQTL_associated_traits.txt') %>% as.data.frame()
color <- fread('Data/color_annotation_v3.txt')
num <- merge(num, color, by.x = 'tissue', by.y = 'tissue_rQTL')
write.table(num, 'Result/03_enrich_GWAS/rhyQTL_associated_traits_fullname.txt', quote = F, row.names = F, sep  = '\t')

theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8))

p1 <- ggplot(data = filter(num, trait.num >= 150), aes(x = trait.num, y = reorder(ShortName, trait.num))) + geom_col(fill = 'red', width = 0.3) + theme +
  labs(x = '# of traits/diseases', y = '')
p1

library(grid)
pdf('./Figure/03_enrich_GWAS/Figure_overlap_with_unique_trait_gt150.pdf', width = 2.3, height = 2.6)
print(p1, vp=viewport(.8, .9, x = 0.5, y = 0.5))
dev.off()



############################################
#
# QQ plot
#
############################################
# Observed
rqtl <- readRDS('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/00_rQTL_mapping/04_Results/cis_rQTL.rds')
GWAS <- fread('./03_enrich_GWAS/GWAS_enrichment/GWAS_r2024-01-19.filter.txt')
#random.all.rqtl <- readRDS('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/02_Functional_annotation_enrichment/all.rQTL.random_sampled_ID.rds')
rsid <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Data/03_MAF059.VCF.MAF.rsid.txt')

obs.overlap <- mclapply(rqtl, function(tmp){
  
  qtl_id <- rsid$rsid[rsid$ID %in% as.character(tmp$ID)]
  gwas_overlap <- GWAS[GWAS$SNP %in% unique(qtl_id),c(1,9)]
  
  return(gwas_overlap)
  
}, mc.cores = 45)

# count the number of overlap
obs.count <- mclapply(obs.overlap, function(tmp){
  
  return(length(unique(tmp$SNP)))
  
}, mc.cores = 45)

obs.count <- as.data.frame(do.call(rbind, obs.count)) %>% rownames_to_column(., var = 'tissue')
names(obs.count) <- c('tissue', 'rqtl.obs')


pval.list <- lapply(obs.overlap, function(tmp){
  
  num <- length(tmp$SNP)
  random_row <- sample(1:nrow(GWAS), num)
  random <- GWAS[random_row,]
  
  pval <- data.frame(
    expect = random[order(random$P_value),]$P_value,
    obs = tmp[order(tmp$P_value),]$P_value) 
  
  return(pval)
})

pval.list <- as.data.frame(do.call(rbind, pval.list)) %>% rownames_to_column(., var = 'tissue')
pval.list$tissue <- gsub('\\.[0-9]+', '', pval.list$tissue)
write.table(pval.list, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/pvalue_correlation/p_value.txt', quote = F, sep = '\t', row.names = F)



############################
#
# plot
#
############################
## get # of rhyQTl
num.list <- lapply(rqtl, function(tmp){
  
  return(length(unique(tmp$ID)))
  
})
names(num.list) <- names(rqtl)
num <- as.data.frame(do.call(rbind, num.list)) %>% rownames_to_column(., var = 'tissue')
names(num) <- c('tissue', 'rhyQTL.num')


theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        #axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8))+
  theme(legend.position = "none")

pval.list <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/pvalue_correlation/p_value.txt')
color <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Data/color_annotation_v3.txt')
num <- num[order(-num$rhyQTL.num),]

############ Top 10
top10 <- filter(pval.list, tissue %in% num$tissue[1:10])
top10$tissue <- factor(top10$tissue, levels = unique(top10$tissue))
p1 <- ggplot(data = top10, 
             aes(x = -log10(expect + 1e-320), y = -log10(obs + 1e-320), color = tissue)) + geom_point(size = 0.1, alpha = 0.8) + 
  geom_abline(intercept = 0, slope = 1) +
  scale_color_manual(values = color$ColorPlot[match(levels(top10$tissue), color$tissue_rQTL)]) +
  theme + 
  #ggtitle('The top 10 tissues with the highest # of rQTLs')+
  #theme(plot.title = element_text(size = 9)) +
  labs(x = 'Expected -log10(P value)', y = 'Obesrved -log10(P value)')
p1

pdf('./Figure/03_enrich_GWAS/p_value_top10_tissue.pdf', height = 1.8, width = 1.8, useDingbats = F)
print(p1)
dev.off()
