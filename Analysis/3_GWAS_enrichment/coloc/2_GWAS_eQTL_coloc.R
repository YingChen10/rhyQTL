
suppressMessages(library(tidyverse))
#suppressMessages(library(stringr))
suppressMessages(library(ieugwasr))
suppressMessages(library(coloc))
suppressMessages(library(data.table))
suppressMessages(library(R.utils))
suppressMessages(library(parallel))


########## FUNCTION: for eQTL only input interested genes ##########
filter_and_read_eQTL <- function(eQTL_file, genes) {
  
  # generate random file name
  temp_file <- tempfile(fileext = ".txt", tmpdir = Temp.dir)
  temp_file2 <- tempfile(fileext = ".txt", tmpdir = Temp.dir)
  genes = c("gene_id", genes)
  
  # write gene name to the temp file
  write.table(genes, temp_file, quote = FALSE, col.names = FALSE, row.names = FALSE)
  
  # grep results of the gene and write to the temp file
  cmd <- paste0("zcat ", eQTL_file, " | grep -F -w -f ", temp_file, " > ", temp_file2)
  system(cmd)
  
  # read results of the gene
  eQTL_filtered <- fread(temp_file2, header = TRUE)
  
  # delete temp files
  unlink(temp_file)
  unlink(temp_file2)
  
  return(eQTL_filtered)
}
####################################################################

######### FUNCTION: format eQTL for coloc ############
Format.eQTL <- function(eqtl_list, maf, genes){
  eqtl_list <- mutate(eqtl_list, Sample.size = round(ma_count/maf)/2)
  eqtl_list <- merge(select(maf, ID, snp, chr, pos, REF, ALT), eqtl_list, by.x = 'ID', by.y = 'variant_id')
  eqtl_list=as.data.frame(eqtl_list)
  eqtl_list$chr <- as.numeric(gsub("chr","",eqtl_list$chr))
  
  eqtl <- filter(eqtl_list,
                 gene_id == genes &
                   chr == tss[tss$Gene == genes,]$Chr &
                   pos >= tss[tss$Gene == genes,]$TSS-1000000 &
                   pos <= tss[tss$Gene == genes,]$TSS+1000000) %>% unique()
  
  eqtl <- transmute(eqtl,
                    ID,
                    gene = gene_id,
                    Sample.size,
                    P.value = pval_nominal,
                    snp, chr, pos, REF, ALT,
                    MAF.Freq = maf, 
                    Beta = slope,
                    SE = slope_se)
  
  eqtl <- filter(eqtl, MAF.Freq > 0) 
  
  return(eqtl)
}
####################################################################

##### FUNCTION: coloc of eQTL and GWAS  ############################
Run.eQTL.GWAS.coloc <- function(eqtl, gwas){
  
  
  if(nrow(eqtl) < 2){return(NULL)}
  
  input2 <- inner_join(eqtl, gwas, by = c("snp","chr", "REF","ALT")) %>% unique()
  
  if(nrow(input2) < 2){return(NULL)}
  
  input2 <- input2[!duplicated(input2$snp),]
  
  if(input2$type[1]=="cc"){
    dataset_gwas2=list(pvalues=input2$pvalues,
                       snp=input2$snp, 
                       type=input2$type[1], 
                       s=input2$s[1],
                       beta=input2$beta,
                       varbeta=input2$varbeta,
                       MAF=input2$MAF,
                       N=input2$N,
                       z=input2$z)
    
  }else{
    dataset_gwas2=list(pvalues=input2$pvalues,
                       snp=input2$snp, 
                       type=input2$type[1], 
                       #s=input$s[1],
                       beta=input2$beta,
                       varbeta=input2$varbeta,
                       MAF=input2$MAF,
                       N=input2$N,
                       z=input2$z
    )
  }
  
  #input maf if GWAS don't have beta and se
  dataset_eQTL <- list(pvalues=input2$P.value, 
                       type="quant", 
                       N=input2$Sample.size,
                       snp=input2$snp,
                       MAF=input2$MAF.Freq,
                       beta=input2$Beta,
                       varbeta=input2$SE)
  
  
  check_dataset(dataset_gwas2)
  check_dataset(dataset_eQTL)
  
  result2 <- coloc.abf(dataset_eQTL,dataset_gwas2)
  
  o2 = data.frame(
    gene = eqtl$gene[1],
    # GWAS=input2$id.y[1],
    eQTL_GWAS.PP0 = result2$summary[2],
    eQTL_GWAS.PP1 = result2$summary[3],
    eQTL_GWAS.PP2 = result2$summary[4],
    eQTL_GWAS.PP3 = result2$summary[5],
    eQTL_GWAS.PP4 = result2$summary[6],
    eQTL_GWAS.SNP_size = nrow(input2))
  
  return(o2)  
} 
####################################################################

dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/"
cat <- "cis_QTL/"
Temp.dir <- paste0(dir, cat, '03_enrich_GWAS/coloc/Temp/')

# creat Temp dir
if (!file.exists(Temp.dir)){dir.create(Temp.dir)}

indir <- '/mnt/disk5_7T/GWASCatalog_summary2/'
outdir <- paste0(dir, cat, 'Result/03_enrich_GWAS/coloc/')

args <- commandArgs(trailingOnly = T)
gwas_id <- args[1]
tissue <- 'Liver'

# file names 
gwas_file <- paste0(dir, cat, "Result/03_enrich_GWAS/coloc/00_Data_GWAS_input/", gwas_id, ".sort.txt.gz")
eQTL_file <- paste0(dir, cat, "03_enrich_GWAS/gwas/gtex.v8.eQTL/", tissue, ".allpairs.txt.gz")

# read GWAS 
gwas <- fread(gwas_file, head=T)
gwas <- as.data.frame(gwas)
gwas$chr <- as.numeric(gwas$chr)
gwas$N <- as.numeric(gwas$N)

# remove duplicated snp
duplicated_rows <- duplicated(gwas$snp)
gwas <- gwas[!duplicated_rows, ]
if(is.na(match('s', names(gwas)))){
  gwas$s <- 0.33
  message('no s')
}
gwas <- gwas[complete.cases(gwas[, c("snp", 'pvalues', "type", 'beta', 's','varbeta', 'MAF', 'N', 'z')]), ]

# read eQTL
egene.list <- fread(paste0(dir, 'eQTLs/GTEx_Analysis_v8_eGenes.txt'))
gene <- unique(egene.list$gene_id[egene.list$Tissue == tissue])
# rqtl_col <- fread(paste0(outdir, '04_rhyQTL_PP4_filter_plot/Liver_04_PP4_filter_plot.txt'))
# gene <- filter(rqtl_col, gwas == gwas_id)$gene

if (length(gene) == 0) {

  message(paste0("no for ", gwas_id, " !"))
  quit(save = "no")

  }


# read tss site 
tss <- read.table(paste0(dir, cat, "03_enrich_GWAS/coloc/Gene.TSS.GRCh38.txt"), head=T, sep="\t")
tss$Chr <- as.numeric(tss$Chr)

# read maf and rsid information
maf <-fread(paste0(dir, 'cis_QTL/Data/03_MAF059.VCF.MAF.rsid.txt'))
names(maf)[names(maf)=='rsid'] <- 'snp'
####################################################################

eqtl.gwas.coloc.result <- lapply(gene, function(x){
  
  # get gene name to test
  genes <- x
  
  # read and filter eQTL file and only keep information of the gene
  eqtl_list <- filter_and_read_eQTL(eQTL_file, genes)
  
  # get eQTL withing TSS flanking 1 Mb region and add rsid information
  eqtl <- Format.eQTL(eqtl_list, maf, genes)
  # coloc of eQTL and GWAS
  eqtl.gwas.coloc <- Run.eQTL.GWAS.coloc(eqtl, gwas)
  
  eqtl.gwas.coloc$GWAS <- gwas_id
  #eqtl.gwas.coloc$rQTL_file <- split_pos_file
  eqtl.gwas.coloc$tissue <- tissue
  
  return(eqtl.gwas.coloc)
  
})

eqtl.gwas.coloc.result <- as.data.frame(do.call(rbind, eqtl.gwas.coloc.result))

outdir = paste0(dir, cat, 'Result/03_enrich_GWAS/coloc/02_eQTL_coloc/')
if (!file.exists(outdir)){dir.create(outdir)}

write.table(eqtl.gwas.coloc.result, paste0(outdir, tissue, '.', gwas_id, '.all.txt'), sep = '\t', quote = F, row.names = F)

