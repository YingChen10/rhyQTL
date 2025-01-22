suppressMessages(library(tidyverse))
#suppressMessages(library(stringr))
#suppressMessages(library(ieugwasr))
suppressMessages(library(coloc))
suppressMessages(library(data.table))
#suppressMessages(library(R.utils))
suppressMessages(library(parallel))
#library(gwasvcf)


###################################### funcation

gwas.format <- function(outdir, gwas_id, indir){
  
  # get  gwas summary harmonised file name
  harmonised <- fread(paste0(dir, cat, '03_enrich_GWAS/GWASCatalog-harmonised-summary-statistics-20240213.txt'), header = F, sep = "\t", na.strings = c("", "NA", "N/A"), quote="")
  names(harmonised) <- 'utrl'
  prefix <- 'https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/'
  harmonised$id <- str_split_fixed(harmonised$utrl, '/', 4)[ , 3]
  harmonised$EFO <- str_split_fixed(harmonised$utrl, '-', 4)[ , 4] %>% gsub('.h.tsv.gz', '', .)
  gwas.filter <- fread('03_enrich_GWAS/GWAS_enrichment/GWAS_r2024-01-19.filter.txt')
  harmonised.filter <- filter(harmonised, id %in% unique(gwas.filter$Study_accession))
  harmonised.filter$file <- str_split_fixed(harmonised.filter$utrl, '/', 5)[ , 5]
  write.table(harmonised.filter, paste0(outdir, 'harmonised-summary.txt'), quote = F, sep = '\t', row.names = F)

  ###################### get sample size
  sample <- fread(paste0(dir, cat, './03_enrich_GWAS/GWAS_enrichment/GWAS_r2024-01-19.filter.txt'))
  sample <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/03_enrich_GWAS/gwas_catalog_v1.0.2.1-studies_r2024-01-19.tsv', sep = "\t", na.strings = c("", "NA", "N/A"), quote="")
  sample <- sample[!is.na(sample$`SUMMARY STATS LOCATION`), ]

  sample$max_value <- sapply(str_extract_all(sample$`INITIAL SAMPLE SIZE`, "\\d+(,\\d+)*"), function(x) {
    # remove ',' and as number
    x_numeric <- as.numeric(gsub(",", "", x))
    # get max value
    max(x_numeric, na.rm = TRUE)
  })

  size <- select(sample, `INITIAL SAMPLE SIZE`, max_value, `STUDY ACCESSION`)
  names(size) <- c('initial.sample.size', 'sample.size', 'id')
  merge <- merge(harmonised.filter, size, by = 'id')

  write.table(merge, paste0(outdir, 'harmonised-summary.txt'), quote = F, sep = '\t', row.names = F)

  harmonised.filter <- fread(paste0(outdir, 'harmonised-summary.txt'))
  ####### impore gwas file
  file <- harmonised.filter$file[harmonised.filter$id == gwas_id]
  gwas <- fread(paste0(indir, file))
  
  
  #### format GWAS results
  if(sum(grepl('hm_rsid', names(gwas))) == 0){
    gwasResults <- select(gwas, rsid, chromosome, base_pair_location, p_value, effect_allele_frequency, beta, standard_error, other_allele, effect_allele)
  }else{
    # if(sum(grepl('hm_beta', names(gwas))) == 0)
    # {
    #   gwasResults <- select(gwas, hm_rsid, hm_chrom, hm_pos, p_value, hm_effect_allele_frequency, beta, standard_error, hm_other_allele, hm_effect_allele)
    # }else{
      gwasResults <- select(gwas, hm_rsid, hm_chrom, hm_pos, p_value, hm_effect_allele_frequency, hm_beta, standard_error, hm_other_allele, hm_effect_allele)
   # }
    
  } 
  
  names(gwasResults) <- c('SNP', 'CHR', 'BP', 'P', 'AF', 'beta', 'SE', 'REF', 'ALT')
  gwasResults$AF[is.na(gwasResults$AF)] <- 0.5
  gwasResults$AF[gwasResults$AF > 0.5] <- 1 - gwasResults$AF[gwasResults$AF > 0.5]
  
  
  S = 0.33
  sample.size <- harmonised.filter$sample.size[harmonised.filter$id == gwas_id]
  type1 <- 'quant'
  
  gwas.table <- gwasResults %>% 
    {list(pvalues = .$P, 
          N = sample.size, 
          MAF = .$AF, 
          beta = .$beta, 
          varbeta = .$SE^2, 
          type = type1, 
          snp = .$SNP, 
          z = .$beta / .$SE, 
          chr = .$CHR, 
          pos = .$BP, 
          REF=.$REF,
          ALT=.$ALT,
          id = gwas_id)
    } %>%  as.data.frame()
  gwas.table$s <- S
  
  gwas.outdir <- paste0(outdir, "/00_Data_GWAS_input/")
  if (!file.exists(gwas.outdir)) {dir.create(gwas.outdir)}
  fwrite(gwas.table, paste0(gwas.outdir, gwas_id, ".sort.txt"), col.names =T , row.names = F, sep="\t", quote = F, na = "NA")
  system(paste0('gzip ',paste0(gwas.outdir, gwas_id, ".sort.txt")))
}
###############################################################





################################################################################### step 1. GWAS input #######################################################################
dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/"
cat <- "cis_QTL/"

indir <- '/mnt/disk5_7T/GWASCatalog_summary2/'
outdir <- paste0(dir, cat, 'Result/03_enrich_GWAS/coloc/')

args <- commandArgs(trailingOnly = T)
gwas_id <- args[1]

file_path <- paste0(outdir, "/00_Data_GWAS_input/", gwas_id, ".sort.txt.gz")

if (file.exists(file_path)) {
  stop(paste0("exits ", gwas_id, " and end!"))
} else {
  gwas.format(outdir, gwas_id, indir)
}

