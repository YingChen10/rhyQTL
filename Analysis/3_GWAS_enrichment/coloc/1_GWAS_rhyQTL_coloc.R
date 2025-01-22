rm(list=ls())
suppressMessages(library(tidyverse))
suppressMessages(library(ieugwasr))
suppressMessages(library(coloc))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))

################################################## function ################################
Run.COLOC <- function(eqtl, gwas, gene_name){
  
  result_out <- data.frame() 
  result_detail <- data.frame()
  result_sig <- data.frame()
  
  eqtl = arrange(eqtl, pos)
  #input = inner_join(eqtl, gwas, by=c("snp","chr","REF","ALT"))
  input = merge(eqtl, gwas, by=c("snp"))  
  
  if(nrow(input)<2){return(NULL)}
  if(input$type[1]=="cc"){
    dataset_gwas=list(pvalues=input$pvalues,
                      snp=input$snp, 
                      type=input$type[1], 
                      s=input$s[1],
                      beta=input$beta,
                      varbeta=input$varbeta,
                      MAF=input$MAF,
                      N=input$N,
                      z=input$z)
    
  }else{
    dataset_gwas=list(pvalues=input$pvalues,
                      snp=input$snp, 
                      type=input$type[1], 
                      #s=input$s[1],
                      beta=input$beta,
                      varbeta=input$varbeta,
                      MAF=input$MAF,
                      N=input$N,
                      z=input$z
    )
  }
  
  # if(min(dataset_gwas$pvalues) > 1e-8){return(NULL)}
  #input maf if GWAS don't have beta and se
  dataset_qtl=list(pvalues=input$p.value, 
                   type="quant", 
                   N=input$Sample.size,
                   snp=input$snp,
                   MAF=input$MAF.Freq)
  
  
  check_dataset(dataset_gwas)
  check_dataset(dataset_qtl)
  
  result <- coloc.abf(dataset_qtl, dataset_gwas)
  o=data.frame(
    gwas.min.p = min(dataset_gwas$pvalues),
    rqtl.min.p = min(dataset_qtl$pvalues),
    rGene=gene_name,
    GWAS=input$id[1],
    PP0=result$summary[2],
    PP1=result$summary[3],
    PP2=result$summary[4],
    PP3=result$summary[5],
    PP4=result$summary[6],
    SNP_size=nrow(input)
  )
  return(o)
}
###########################################################################################

################################################################################### step 2. GWAS rhyQTL coloc #######################################################################
dir <- "/N/scratch/cc123/gwas/rQTL/"
#dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/"
cat <- "cis_QTL/"
split_pos_file <- list.files(paste0(dir, cat, '00_rQTL_mapping/04_hanova/Liver/')) %>% gsub('.rds', '', .)


args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
gwas_id <- args[2]

message(paste0('processing ', gwas_id, '!'))

outdir <- paste0(dir, cat, 'Result/03_enrich_GWAS/coloc/01_rhyQTL_coloc/', tissue)
outdir.coloc <- paste0(dir, cat, 'Result/03_enrich_GWAS/coloc/01_rhyQTL_coloc_PP4_filter/', tissue)


if (!file.exists(outdir)) {dir.create(outdir, recursive = T)}
if (!file.exists(paste0(outdir, gwas_id))) {dir.create(paste0(outdir, gwas_id), recursive = T)}
if (!file.exists(outdir.coloc)) {dir.create(outdir.coloc, recursive = T)}


if (file.exists(paste0(outdir, gwas_id, '.txt'))) {
message(paste0("exits ", gwas_id, " and end!"))
quit(save = "no")
}

## read GWAS #################################################
gwas_file = paste0(dir, cat, "Result/03_enrich_GWAS/coloc/00_Data_GWAS_input/", gwas_id, ".sort.txt.gz")
gwas = fread(gwas_file, head = T)
gwas$chr <- as.numeric(gwas$chr)
gwas$pvalues <- as.numeric(gwas$pvalues)
# remove duplicated snp
duplicated_rows <- duplicated(gwas$snp)
gwas <- gwas[!duplicated_rows, ]
if(is.na(match('s', names(gwas)))){
  gwas$s <- 0.33
  message('no s')
  }
gwas <- gwas[complete.cases(gwas[, c("snp", 'pvalues', "type", 'beta', 's','varbeta', 'MAF', 'N', 'z')]), ]

# read maf and rsid information
maf <-fread(paste0(dir, cat, 'Data/03_MAF059.VCF.MAF.rsid.txt'))
names(maf)[names(maf)=='rsid'] <- 'snp'

all <- lapply(split_pos_file, function(split_pos){

  if (file.exists(paste0(outdir, gwas_id, '/',split_pos, '.txt'))) {
    message(paste0("exits ", gwas_id, ' ',split_pos, " and end!"))
    quit(save = "no")
  }
  
  rqtl_file <- paste0(dir, cat, '00_rQTL_mapping/04_hanova/', tissue, '/', split_pos, ".rds")
  rqtl_list <- readRDS(rqtl_file)
  
  regression_file <- paste0(dir, cat, '00_rQTL_mapping/01_Rhythm_regression/', tissue, '/', split_pos, ".rds")
  regression.tmp <- readRDS(regression_file)
  
  intersect <- intersect(names(rqtl_list), names(regression.tmp))
  
  split_file_coloc <- lapply(1:length(intersect), function(i) {
    #split_file_coloc <- lapply(1:2, function(i) {  
    
    print(i)
    gene_name <- intersect[i]  
    message(gene_name)
    gene.tmp <- rqtl_list[[gene_name]]           
    
    #names(gene.tmp) <- c('ID', 'geno.0','geno.1','geno.2','geno.0.sample','geno.1.sample','geno.2.sample','HANOVA.no.norm',  'HANOVA.norm')
    gene.tmp$HANOVA.norm <- as.numeric(gene.tmp$HANOVA.norm)
    
    rqtl <- merge(gene.tmp, maf, by = 'ID')  
    duplicated_rows <- duplicated(rqtl$snp)
    rqtl <- rqtl[!duplicated_rows, ]
    
    rqtl$Sample.size <- as.numeric(rqtl$geno.0.num) + as.numeric(rqtl$geno.0.num)  
    names(rqtl)[names(rqtl) == 'HANOVA.norm'] <- 'p.value'
    rqtl$p.value[rqtl$p.value == 0] <- 1e-20
    coloc <- Run.COLOC(rqtl, gwas, gene_name)
    message(coloc)
    return(coloc)
  })
  
  split_file_coloc <- as.data.frame(do.call(rbind, split_file_coloc))
  split_file_coloc$file <- split_pos
  
  write.table(split_file_coloc, paste0(outdir, gwas_id, '/',  split_pos, '.txt'), sep = '\t', quote = F, row.names = F)
  
  
  return(split_file_coloc)

})


all <- as.data.frame(do.call(rbind, all))
write.table(all, paste0(outdir, gwas_id, '.txt'), sep = '\t', quote = F, row.names = F)

coloc_sig <- filter(all, PP4 >= 0.5)
if(nrow(coloc_sig)>0){write.table(coloc_sig, paste0(outdir.coloc, gwas_id, '.sig.txt'), sep = '\t', quote = F, row.names = F)}

