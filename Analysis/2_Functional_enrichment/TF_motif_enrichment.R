library(data.table)
library(tidyverse)

setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')
##############################################

# generate motif and GTEx snp relationship

##############################################
db <- read.table('02_Functional_annotation_enrichment/motif_enrichment_OR/homer_motif_db.txt')
names(db) <- c('motif', 'name')
db$TF <- gsub('/.*/Homer', '', db$name)

# excluse redundant motifs
uni.motif <- db[!duplicated(db$motif),]
uni.motif$TF <- gsub("\\([^\\)]+\\)", "", uni.motif$TF)
uni.motif <- uni.motif[!duplicated(uni.motif$TF),]
db <- uni.motif
write.table(db, '02_Functional_annotation_enrichment/motif_enrichment_OR/homer_motif_db_mod.txt', quote = F, sep = '\t', row.names = F)

# excluse redundant TFs
uni.motif <- filter(uni.motif, name != 'c-Myc(bHLH)/LNCAP-cMyc-ChIP-Seq(Unpublished)/Homer' &
         name != 'COUP-TFII(NR)/K562-NR2F1-ChIP-Seq(Encode)/Homer' &
         name != 'GRE(NR),IR3/RAW264.7-GRE-ChIP-Seq(Unpublished)/Homer' &
         name != 'RE(HSF)/Striatum-HSF1-ChIP-Seq(GSE38000)/Homer' &
         name != 'Nr5a2(NR)/mES-Nr5a2-ChIP-Seq(GSE19019)/Homer' &
         name != 'OCT:OCT(POU,Homeobox)/NPC-Brn1-ChIP-Seq(GSE35496)/Homer' &
         name != 'STAT6(Stat)/Macrophage-Stat6-ChIP-Seq(GSE38377)/Homer' &
         name != 'p53(p53)/mES-cMyc-ChIP-Seq(GSE11431)/Homer')

unique(uni.motif$TF[duplicated(uni.motif$TF)])

uni.motif <- filter(uni.motif,
                    name != 'HRE(HSF)/Striatum-HSF1-ChIP-Seq(GSE38000)/Homer' &
                      name != 'THRb(NR)/Liver-NR1A2-ChIP-Seq(GSE52613)/Homer')



motif_input <- '02_Functional_annotation_enrichment/motif_enrichment_OR/SNP_motif_homer/chr'
motif <- list()
i <- 1
motif.list <- list()
motif <- fread(paste0(motif_input, i, '_snp.bed')) %>% as.data.frame()
for(x in unique(db$motif)){

  motif.list[[x]] <- filter(motif, V8 == x)[,c('V8', 'V4')]

  }


for(i in 2:22){

  motif <- fread(paste0(motif_input, i, '_snp.bed')) %>% as.data.frame()

  for(x in unique(db$motif)){

    motif.list[[x]] <- rbind(motif.list[[x]], filter(motif, V8 == x)[,c('V8', 'V4')])

  }

}
saveRDS(motif.list, './02_Functional_annotation_enrichment/motif_enrichment_OR/homer.motif.gtex.snp.rds')


##############################################

# calculate enrichment

##############################################
motif.list <- readRDS('./02_Functional_annotation_enrichment/motif_enrichment_OR/homer.motif.gtex.snp.rds')
rqtl <- readRDS('00_rQTL_mapping/04_Results/cis_rQTL.rds')

# get the # of rQTL and motif overlap 
library(parallel)
result <- lapply(rqtl, function(tmp){
  
  or.tissue <- lapply(motif.list, function(motif.tmp) {
    
    if(nrow(motif.tmp) >= 1){
      rqtl.num <- length(unique(tmp$ID))
      motif.rqtl.num <- length(intersect(unique(motif.tmp$V4), unique(tmp$ID)))
      motif.num <- length(unique(motif.tmp$V4))
      
      or <- data.frame(
        motif = unique(motif.tmp$V8),
        rqtl.num = rqtl.num,
        motif.rqtl.num = motif.rqtl.num,
        motif.num = motif.num)
      
      return(or)}
    
  })
  
  return(or.tissue)
  
})
names(result) <- names(rqtl)
saveRDS(result, './Result/02_Functional_annotation_enrichment/motif_enrichment_OR/homer.snp.num.rds')

## calculate odds ratio 
total.num <- 10770861
OR.list <- lapply(result, function(tissue){
  
  
  tmp <- as.data.frame(do.call(rbind, tissue))
  
  OR <- c()
  
  for(i in 1:nrow(tmp)){
    
    num <- c(tmp[i,3], tmp[i,2]-tmp[i,3],  tmp[i,4]-tmp[i,3], 
             total.num - tmp[i,2]-tmp[i,4] + tmp[i,3])
    
    data <- matrix(num, nrow = 2) 
    # colnames(data) <- c("rqtl", "non_rqtl") 
    # rownames(data) <- c("motif", "non_motif")
    
    OR <-data.frame(motif = tmp[i,1], 
                    rqtl.motif = tmp[i,3], 
                    rqtl.non.motif = tmp[i,2]-tmp[i,3],  
                    non.rqtl.motif = tmp[i,4]-tmp[i,3], 
                    non.rqtl.non.motif = total.num - tmp[i,2]- tmp[i,4] + tmp[i,3],
                    odds.ratio = as.numeric(fisher.test(data)$estimate),
                    p.value = as.numeric(fisher.test(data)$p.value)) %>% rbind(OR, .)
    
  }
  
  return(OR)
})

or <- as.data.frame(do.call(rbind, OR.list))
or <- rownames_to_column(or, var = 'tissue')
or$tissue <- gsub('\\.[0-9]+$', '', or$tissue)
db <- fread('02_Functional_annotation_enrichment/motif_enrichment_OR/homer_motif_db_mod.txt')
or <- merge(db, or, by = 'motif')

write.table(or, 'Result/02_Functional_annotation_enrichment/motif_enrichment_OR/homer.or.tissue.txt', quote = F, sep = '\t', row.names = F)



