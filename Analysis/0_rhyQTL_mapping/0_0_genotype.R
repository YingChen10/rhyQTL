rm(list = ls()) 
print (Sys.time())
library(data.table)
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library("dryR"))
suppressMessages(library("lmtest"))

################## parameters ###########################
dir <- "/N/scratch/cc123/gwas/rQTL/"
cat <- "cis_QTL/"

args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
pos_file <- args[2]

# tissue <- 'Brain-Substantianigra'
# pos_file <- 'aa'
outdir <- paste0(dir, cat, '00_rQTL_mapping/00_Genotype/', tissue, '/')
logoutdir <- paste0(dir, cat, 'Log/00_Genotype/', tissue, '/')
if (!file.exists(outdir)) {dir.create(outdir, recursive = TRUE)}
if (!file.exists(logoutdir)) {dir.create(logoutdir, recursive = TRUE)}


a <- Sys.time()
write.table(a, paste0(logoutdir, pos_file), quote = F, row.names = F, col.names = F)

min_size <- 50
period <- 24
group  <- 2 # keep snv with sample size gt  min_size in at least this group number
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
  
  # only keep wt group with > min size 
  num_cp <- filter(num, Freq >= min_size)
  
  if(nrow(num_cp) >= 2){
    # subset geno
    geno <- filter(geno, genotype %in% as.character(num$Var1))
  
    # return 
    result <- list(geno = select(geno, genotype))
    return(result)
  }
}


###### import files ##############
### input time information
time <- fread(paste0(dir, 'GTEx_donor_time_science.txt'))

###### input gene position
pos <- read.table(paste0(dir, cat, "split_pos/", pos_file))
pos$tss <- pos$V2
pos$tss[pos$V4 == '-'] <- pos$V3[pos$V4 == '-']
pos$tss_up <- pos$tss - 1000000
pos$tss_down <- pos$tss + 1000000
pos$tss_up[pos$tss_up < 0] <- 0

### head of snv file
head <- read.table(paste0(dir, 'head_sub.txt'))

# import expression file
expression <- read.table(paste0(dir,'00_data/CPM_covariate_remove/', tissue, '.txt'), header = T)
names(expression) <- str_split_fixed(names(expression), '\\.', 2)[ ,2]
expression$EnsemblID <- str_split_fixed(row.names(expression), '_', 2)[ ,1]


########## retrieve genotype information ##############################
all_results <- list()
for(line in  1:nrow(pos)){
  print (paste0('Starting processing gene', line))
  ## retrieve TSS-1M flanking region  for each gene
  gene_info <- unlist(pos[line, c(1,12,13)])
  region <- paste0(gene_info[1], ':', gene_info[2], '-', gene_info[3])
  gene <- pos$V5[line]
  

  # retrieve SNVs in TSS-1M flanking region  for each gene
  output <- system(paste0("/geode2/home/u020/cc123/BigRed200/miniconda3/bin/tabix ", dir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hwe_MAF059.vcf.gz ", region), intern=TRUE)
  # output <- system(paste0("/workspace/rsrch1/ychen/miniconda3/bin/tabix ", dir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hwe_MAF059.vcf.gz ", region), intern=TRUE)
 
  snv <- read.table(textConnection(output), header=F, sep="\t")
  names(snv) <- as.character(head)
  snv <- snv[grepl("A|G|C|T", snv$REF, perl = T), ]
  snv <- snv[grepl("A|G|C|T", snv$ALT, perl = T), ]
  if(nrow(snv) == 0){
    #   print ('no cis SNVs')
    #all_results[[line]] <- NULL
    next}
  #print('1. Retrieving SNVs in cis region done!')
  
  
  ## retrieve gene expression matrix for the gene
  expr <- filter(expression, EnsemblID == gene) 
  if(nrow(expr) == 0){ #   print ('no gene in expr file')
    #all_results[[line]] <- NULL
    next}
  # remove sample with expression of the gene didn't detected
  expr <- select(expr, names(expr)[!is.na(expr)])
  # print('2. Retrieving expression matrix done!')
  
  
  ## filter snv according samples in expression matrix
  snv_interscet <- select(snv, ID, intersect(names(expr)[-1], names(snv)))
  
  geno_list <- lapply(1:nrow(snv_interscet), Filter_SNV)
  names(geno_list) <- snv_interscet$ID
  geno_list <- geno_list[!sapply(geno_list, is.null)]
  if(length(names(geno_list)) == 0){
    #   print ('no cis SNVs after filter')
    #all_results[[line]] <- NULL
    next}
  
  all_results[[gene]] <- geno_list
}

all_results <-  all_results[!sapply(all_results, is.null)]
saveRDS(all_results, paste0(outdir, pos_file, '.rds'))

b <- data.frame(tissue = tissue, file = pos_file, s = a, e = Sys.time())
write.table(b, paste0(logoutdir, pos_file, '.complete'), sep = '\t', quote = F, row.names = F, col.names = F)