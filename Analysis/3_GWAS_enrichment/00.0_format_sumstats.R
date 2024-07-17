rm(list=ls())
setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/')

suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(tidyverse))

# import which trait ID
# We test in 15 traits ("GCST000755", "GCST002216", "GCST90092808", "GCST90092809", "GCST90092826", "GCST90092875", "GCST90092928", "GCST90092975", "GCST90092980", "GCST90092981", "GCST90092985", "GCST90092987", "GCST90093002", "GCST90269502", "GCST90269582")
args <- commandArgs(trailingOnly = T)
trait.ID <- args[1]
# trait.ID <- 'GCST90018975'


term <- read.table('03_enrich_GWAS/GWAScatalog_201401/gwas.parent.term.list', header = T)

# file path information at local
gwas_indir <- '/mnt/disk/GWASCatalog_summary/'


# # download gwas summary harmonised
harmonised <- fread('03_enrich_GWAS/GWASCatalog-harmonised-summary-statistics-20240213.txt', header = F, sep = "\t", na.strings = c("", "NA", "N/A"), quote="")
names(harmonised) <- 'utrl'
prefix <- 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
harmonised$id <- str_split_fixed(harmonised$utrl, '/', 4)[ , 3]

gwas.filter <- fread('03_enrich_GWAS/GWAS_enrichment/GWAS_r2024-01-19.filter.txt')
gwas.filter.unique <- unique(select(gwas.filter, Study_accession, Trait_mod))

harmonised.filter <- filter(harmonised, id %in% unique(gwas.filter$Study_accession))
harmonised.filter$file_name <- str_split_fixed(harmonised.filter$utrl, '/', 5)[ , 5]
harmonised.filter <- merge(harmonised.filter, gwas.filter.unique, by.x = 'id', by.y = 'Study_accession')


## import sample size 
GWAS.association <- fread('03_enrich_GWAS/gwas_catalog_v1.0.2-associations_e111_r2024-01-19.tsv', header = T, sep = "\t", na.strings = c("", "NA", "N/A"), quote="")
names(GWAS.association) <- gsub(' ', '_', names(GWAS.association))
GWAS.association$sample.size <- as.numeric(gsub(",", "", str_extract(GWAS.association$INITIAL_SAMPLE_SIZE, "\\b[\\d,]+\\b")))
a <- select(GWAS.association, sample.size, INITIAL_SAMPLE_SIZE)
sample.size <- unique(select(GWAS.association, STUDY_ACCESSION, sample.size))


# outdir
outdir <- '/mnt/disk/LDSC/sumstats/'

file.path <- paste0(gwas_indir, harmonised.filter$file_name[harmonised.filter$id == trait.ID])


gwas <- fread(file.path)

if(sum(grepl('hm_rsid', names(gwas))) == 0){
  
  sumstats <- data.frame(SNP = gwas$rsid,
                         N = sample.size$sample.size[sample.size$STUDY_ACCESSION == trait.ID],
                         Z = gwas$beta/gwas$standard_error,
                         P = gwas$p_value,
                         A1 = gwas$effect_allele,
                         A2 = gwas$other_allele)
}else{
  
  sumstats <- data.frame(SNP = gwas$hm_rsid,
                         N = sample.size$sample.size[sample.size$STUDY_ACCESSION == trait.ID],
                         Z = gwas$hm_beta/gwas$standard_error,
                         P = gwas$p_value,
                         A1 = gwas$hm_effect_allele,
                         A2 = gwas$hm_other_allele)
  
} 

sumstats <- sumstats %>% mutate(N = as.numeric(N),
                                Z = as.numeric(Z), 
                                P = as.numeric(P)) %>% na.omit(.) %>%
  filter(., P >= 0 & P <= 1)

write.table(sumstats, paste0(outdir, x, '.', harmonised.filter$Trait_mod[harmonised.filter$id == x], '.sumstats'), sep = '\t', row.names = F, quote = F)