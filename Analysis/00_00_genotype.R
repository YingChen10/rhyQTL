rm(list = ls()) 
print (Sys.time())
library(data.table)
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library("dryR"))
suppressMessages(library("lmtest"))

################## parameters ###########################
dir <- "Projects/rQTL/"
cat <- "cis_QTL/"

args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
pos_file <- args[2]


if (!file.exists(paste0(dir, cat, '00_rQTL_mapping_00_Genotype/', tissue))) {
  dir.create(paste0(dir, cat, '00_rQTL_mapping_00_Genotype/', tissue))}


min_size <- 50
period <- 24
group  <- 2 # keep snp with sample size greater than min_size in at least 2 genotype groups
#######################################################

#################### functions ########################
## filter SNV and generate genotype list and downsampling 
Filter_SNV <- function(j){
  
  # calculate sample size in each genotype group
  geno <- as.data.frame(t(snv_interscet[j, 2:ncol(snv_interscet)]))
  names(geno) <- 'V1'
  geno <- subset(geno, geno$V1 != '.')
  geno$genotype <- as.numeric(substring(geno$V1,1,1)) + as.numeric(substring(geno$V1,3,3))
  num <- as.data.frame(table(geno$genotype))
  
  # only keep genotype groups with sample size >= min_size 
  num_cp <- filter(num, Freq >= min_size)
  
  if(nrow(num_cp) >= 2){
    # subset samples with their genotype
    geno <- filter(geno, genotype %in% as.character(num$Var1))
  
    # return 
    result <- list(geno = select(geno, genotype))
    return(result)
  }
}


###### import files ##############
### time information
time <- fread(paste0(dir, 'GTEx_donor_time.txt'))

### gene position file, TSS information of 18,200 protein coding genes are split into 49 files in order to reduce the calculation time
pos <- read.table(paste0(dir, cat, "split_pos/", pos_file))
pos$tss <- pos$V2
pos$tss[pos$V4 == '-'] <- pos$V3[pos$V4 == '-']
pos$tss_up <- pos$tss - 1000000
pos$tss_down <- pos$tss + 1000000
pos$tss_up[pos$tss_up < 0] <- 0

### head of snv file
head <- read.table(paste0(dir, 'head_sub.txt'))

### input expression file
expression <- fread(paste0(dir,'GTEx_nor_expression/',tissue, '.txt'))
expression <- as.data.frame(expression)
expression <- expression[,1:(ncol(expression) - 34)]
expression$EnsemblID <- sapply(strsplit(expression$EnsemblID, '_'), "[", 1)
names(expression)[2:ncol(expression)] <- sapply(strsplit(names(expression)[2:ncol(expression)], '\\.'), "[", 2)

########## retrieve genotype information at each cis-variant locus for each gene #########################
all_results <- list()
for(line in  1:nrow(pos)){
  print (paste0('Starting processing gene', line))
  
  # get TSS-1M flanking region for the gene
  gene_info <- unlist(pos[line, c(1,12,13)])
  region <- paste0(gene_info[1], ':', gene_info[2], '-', gene_info[3])
  gene <- pos$V5[line]
  

  # retrieve SNVs within TSS-1M flanking region for the gene
  output <- system(paste0("tabix ", dir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hwe_MAF059.vcf.gz ", region), intern=TRUE)
 
  snv <- read.table(textConnection(output), header=F, sep="\t")
  names(snv) <- as.character(head)
  snv <- snv[grepl("A|G|C|T", snv$REF, perl = T), ]
  snv <- snv[grepl("A|G|C|T", snv$ALT, perl = T), ]
  if(nrow(snv) == 0){
    
    print ('no cis SNVs')
    #all_results[[line]] <- NULL
    next}
  
  print('1. Retrieving SNVs in cis region done!')
  
  
  # retrieve gene expression matrix for the gene
  expr <- filter(expression, EnsemblID == gene) 
  if(nrow(expr) == 0){ 
    print ('no gene in expression file')
    #all_results[[line]] <- NULL
    next}
  
  # remove samples in which the gene expression was not detected
  expr <- select(expr, names(expr)[!is.na(expr)])
  print('2. Retrieving expression matrix done!')
  
  
  ## filter snv according to the sample id in expression matrix
  snv_interscet <- select(snv, ID, intersect(names(expr)[-1], names(snv)))
  
  geno_list <- lapply(1:nrow(snv_interscet), Filter_SNV)
  names(geno_list) <- snv_interscet$ID
  geno_list <- geno_list[!sapply(geno_list, is.null)]
  if(length(names(geno_list)) == 0){
    print ('After filtering, there are no cis SNVs remaining')
    #all_results[[line]] <- NULL
    next}
  
  all_results[[gene]] <- geno_list
}

all_results <-  all_results[!sapply(all_results, is.null)]

saveRDS(all_results, paste0(dir, cat, '00_rQTL_mapping_00_Genotype/', tissue, '/', pos_file, '.rds'))
print (Sys.time())