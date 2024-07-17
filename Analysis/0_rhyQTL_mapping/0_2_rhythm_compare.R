############################################### 

#####   2. downsample and compare rhythm  #####

###############################################


rm(list = ls()) 

library(data.table)
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library("dryR"))
suppressMessages(library("lmtest"))

############################# parameters #########################
dir <- "rQTL/"
cat <- "cis_QTL/00_rQTL_mapping/"

args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
pos_file <- args[2]

indir <- paste0(dir, cat, '01_Rhythm_regression_filter/', tissue, '/')
outdir <- paste0(dir, cat, '02_Rhythm_compare/', tissue, '/')
logoutdir <- paste0(dir, 'cis_QTL/Log/02_Rhythm_compare/', tissue, '/')

if (!file.exists(outdir)) {dir.create(outdir)}
if (!file.exists(logoutdir)) {dir.create(logoutdir)}

# tissue <- 'Brain-Hippocampus'
# pos_file <- 'split_pos_aa'

start <- Sys.time()
write.table(start, paste0(logoutdir, pos_file), sep = '\t', quote = F, row.names = F, col.names = F)

############################ function ###########################
# The genotype group with the more samples was subsampled to match the sample size of the smaller group
random.sample.test <- function(tmp, num){
  
  random_sample <- sample(row.names(filter(tmp, genotype == num$Var1[1])), size = num$Freq[2], replace = FALSE)
  random_tmp <- rbind(tmp[random_sample,], filter(tmp, genotype == num$Var1[2]))
  dryList = drylm(random_tmp$x, random_tmp$genotype, random_tmp$time)
  parameters <- dryList[["parameters"]][1,]
  return(parameters)
  
}
##################################################################


#################### import files ################################
# import genotype information
all_results <-  readRDS(paste0(indir, pos_file, '.geno.rds'))

# import thy harmonic regression fit result
filter.result <- readRDS(paste0(indir, pos_file, '.rds'))

# input time information
time <- fread(paste0(dir, 'GTEx_donor_time_science.txt'))

# # input expression file
expression <- fread(paste0(dir,'GTEx_nor_expression/',tissue, '.txt'))
expression <- as.data.frame(expression)
expression <- expression[,1:(ncol(expression) - 34)]
expression$EnsemblID <- sapply(strsplit(expression$EnsemblID, '_'), "[", 1)
names(expression)[2:ncol(expression)] <- sapply(strsplit(names(expression)[2:ncol(expression)], '\\.'), "[", 2)


######################## Main ####################################
dryR.results.list <- list()
for(i in 1:length(filter.result)){

  gene <- names(filter.result)[i]
  fit.filter <- filter.result[[gene]]

  if(nrow(fit.filter) == 0){next}

  dif.rhy.list <- lapply(row.names(fit.filter), function(id){

    dif.rhy <- c()

    # get genotype information for current SNP in the gene
    geno <- all_results[[gene]][[id]]
    num <- as.data.frame(table(factor(geno$genotype, levels = 0:2)))

    dif.rhy <- c(id, num$Freq)
    num <- num[order(-num$Freq), ]
    geno_tmp <- filter(geno, genotype == num$Var1[1] | genotype == num$Var1[2])

    # organize as data frame to keep time and expression information
    tmp <- data.frame(genotype = geno_tmp$genotype,
                      time = time$hour[match(rownames(geno_tmp), time$SUBJ.ID)],
                      x = as.numeric(select(filter(expression, EnsemblID == gene), row.names(geno_tmp))))
    row.names(tmp) <- row.names(geno_tmp)
    tmp$genotype <- as.character(tmp$genotype)

    # randomly sample to do DryR test
    dif.rhy <- c(dif.rhy, as.numeric(random.sample.test(tmp, num)))

    return(dif.rhy)
  })

  dif.rhy.list <- as.data.frame(do.call(rbind, dif.rhy.list))
  names(dif.rhy.list) <- c('ID', 'sample.size.0', 'sample.size.1', 'sample.size.2', paste0(c('mean', 'a', 'b', 'amp', 'relamp', 'phase'), '_geno1'),
                           paste0(c('mean', 'a', 'b', 'amp', 'relamp', 'phase'), '_geno2'), 'chosen_model', 'chosen_model_BICW', 'chosen_model_mean', 'chosen_model_mean_BICW')

  dif.rhy.list$Gene <- gene
  dif.rhy.list <- dif.rhy.list %>% relocate(Gene, .before = 1)

  dryR.results.list[[gene]] <- dif.rhy.list
}

dryR.results <- as.data.frame(do.call(rbind, dryR.results.list))
dryR.results <- filter(dryR.results, chosen_model != 1 & chosen_model != 4)
write.table(dryR.results, paste0(outdir, pos_file), sep = '\t', quote = F, row.names = F)

a <- data.frame(tissue = tissue, file = pos_file, s = start, n = length(all_results), e = Sys.time())

write.table(a, paste0(logoutdir, pos_file), sep = '\t', quote = F, row.names = F, col.names = F)



