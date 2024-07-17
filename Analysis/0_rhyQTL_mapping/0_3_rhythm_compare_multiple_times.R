
############################################################ 

#####   3. downsample mutiple times to compare rhythm  #####

############################################################


rm(list = ls()) 
library(data.table)
suppressMessages(library(tidyverse))
#suppressMessages(library(dplyr))
suppressMessages(library("dryR"))
#suppressMessages(library("lmtest"))
library(DescTools)

############################# parameters #########################
dir <- "Projects/rQTL/"
cat <- "cis_QTL/00_rQTL_mapping/"

args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
pos_file <- args[2]

# e.g.
# tissue <- 'Brain-Hippocampus'
# pos_file <- 'split_pos_aa'

indir1 <- paste0(dir, cat, '01_Rhythm_regression_filter/', tissue, '/')
indir2 <- paste0(dir, cat, '02_Rhythm_compare/', tissue, '/')
outdir <- paste0(dir, cat, '03_Rhythm_compare_multiple_times', '/', tissue, '/')
logoutdir <- paste0(dir, 'cis_QTL/Log/03_Rhythm_compare_multiple_times/', tissue, '/')

if (!file.exists(paste0(dir, cat, '03_Rhythm_compare_multiple_times/'))) {dir.create(paste0(dir, cat, '03_Rhythm_compare_multiple_times/'))}
if (!file.exists(outdir)) {dir.create(outdir)}
if (!file.exists(paste0(dir, 'cis_QTL/Log/03_Rhythm_compare_multiple_times/'))) {dir.create(paste0(dir, 'cis_QTL/Log/03_Rhythm_compare_multiple_times/'))}
if (!file.exists(logoutdir)) {dir.create(logoutdir)}

start <- Sys.time()
write.table(start, paste0(logoutdir, pos_file), sep = '\t', quote = F, row.names = F, col.names = F)




############################ function ###########################
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
all_results <-  readRDS(paste0(indir1, pos_file, '.geno.rds'))

# import thy fit filter
filter.result <- read.table(paste0(indir2, pos_file), header = T)

# input time information
time <- fread(paste0(dir, 'GTEx_donor_time_science.txt'))

# # input expression file
expression <- fread(paste0(dir,'GTEx_nor_expression/',tissue, '.txt'))
expression <- as.data.frame(expression)
expression <- expression[,1:(ncol(expression) - 34)]
expression$EnsemblID <- sapply(strsplit(expression$EnsemblID, '_'), "[", 1)
names(expression)[2:ncol(expression)] <- sapply(strsplit(names(expression)[2:ncol(expression)], '\\.'), "[", 2)


######################## Main ###############################

compare.list <- lapply(1:nrow(filter.result), function(i){

  gene <- filter.result$Gene[i]
  id <- filter.result$ID[i]
  
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
  random.list <- list()
  
  for(j in 1:20){random.list[[j]] <- as.numeric(random.sample.test(tmp, num))}
  
  random <- as.data.frame(do.call(rbind, random.list))
  names(random) <- c(paste0(c('mean', 'a', 'b', 'amp', 'relamp', 'phase'), '_geno1'),
                     paste0(c('mean', 'a', 'b', 'amp', 'relamp', 'phase'), '_geno2'), 
                     'chosen_model', 'chosen_model_BICW', 'chosen_model_mean', 'chosen_model_mean_BICW')
    
  # G test to get p value and chosen_model
  observed <- as.data.frame(table(random$chosen_model))
  if(nrow(observed) == 1){P = 0}else{
    P <- GTest(observed$Freq, p = rep(1/nrow(observed), nrow(observed)))$p.value}
  
  model.result <- c()
  if(P < 0.05){
    
    model <- observed$Var1[which.max(observed$Freq)]
    
    model.result <- data.frame(
      gene = gene,
      ID = id,
      chosen_model = model,
      gtest.p.value = P,
      time = max(observed$Freq),
      amp_geno1.mean = mean(filter(random, chosen_model == model)$amp_geno1),
      amp_geno2.mean = mean(filter(random, chosen_model == model)$amp_geno2))
    }

  return(model.result)

  })

   
compare.result <- as.data.frame(do.call(rbind, compare.list))
write.table(compare.result, paste0(outdir, pos_file), sep = '\t', quote = F, row.names = F)

# exclude chosen_model 1 and 4 
write.table(filter(compare.result, chosen_model != 1 & chosen_model != 4), paste0(outdir, pos_file, '.filter'), sep = '\t', quote = F, row.names = F)


a <- data.frame(tissue = tissue, file = pos_file, s = start, n = nrow(filter.result), e = Sys.time())
write.table(a, paste0(logoutdir, pos_file), sep = '\t', quote = F, row.names = F, col.names = F)
