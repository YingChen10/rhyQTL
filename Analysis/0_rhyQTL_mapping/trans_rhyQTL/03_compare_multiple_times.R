print (Sys.time())

library(data.table)
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library("dryR"))
suppressMessages(library("lmtest"))

################## parameters ###########################
dir <- '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/'
cat <- "cis_QTL/"

args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
clock_gene <- args[2]

# tissue <- 'Brain-Amygdala'
# clock_gene <- 'NR1D1'

# create out dir
outdir <- paste0(dir, cat, '13_Results/08_trans_rhyQTL/', tissue, '/')
logoutdir <- paste0(dir, cat, 'Log/08_trans_rhyQTL/', tissue, '/')
if (!file.exists(outdir)) {dir.create(outdir, recursive = T)}
if (!file.exists(logoutdir)) {dir.create(logoutdir, recursive = T)}

start <- Sys.time()
write.table(start, paste0(logoutdir, tissue), sep = '\t', quote = F, row.names = F, col.names = F)

min_size <- 50
period <- 24
group  <- 2

############################### FUNCTION #####################################
## downsample
DownSample <- function(geno){
  size2test <- min(as.numeric(table(geno$genotype)))
  idx <- c()
  for(x in unique(geno$genotype)){
    idx_tmp <- sample(row.names(geno[geno$genotype == x,]), size2test)
    idx <- c(idx, idx_tmp)}
  geno = geno[idx, ]
  
  return(geno)
  
}


############################## Main #######################################
### input time information
time <- fread(paste0(dir, 'GTEx_donor_time_science.txt'))

### head of snv file
head <- read.table(paste0(dir, 'head_sub.txt'))

###### input genotype file
info <- readRDS(paste0(outdir, clock_gene, '.01.sample.rds'))
info <- info[!sapply(info, is.null)]
compare <- readRDS(paste0(outdir, clock_gene, '.02.compare.rds'))

# input expression file
expression <- fread(paste0(dir,'GTEx_nor_expression/', tissue, '.txt'))
expression <- as.data.frame(expression)
expression <- expression[,1:(ncol(expression) - 34)]
expression$EnsemblID <- sapply(strsplit(expression$EnsemblID, '_'), "[", 1)
names(expression)[2:ncol(expression)] <- sapply(strsplit(names(expression)[2:ncol(expression)], '\\.'), "[", 2)
expression <- column_to_rownames(expression, var = 'EnsemblID')


print('rQTL mapping statring!')

r <- lapply(1:length(compare), function(i){
  
  x <- info[[i]]
  
  ID <- names(compare)[i]
  
  ### get significant results and Run dryR
  dryR_test <- unique(compare[[i]]$gene)
  
  
compare <- data.frame()  
for(j in 1:20){
   
   # random select
   geno <- DownSample(x$geno)
   expression.test <- na.omit(expression[dryR_test,])
   dryList = drylm(select(expression.test, row.names(geno)), as.character(geno$genotype), time$hour[match(row.names(geno), time$SUBJ.ID)])
   
   dryResults <- dryList[["parameters"]]
   
   compare <- data.frame(gene = row.names(dryResults),
                         chosen_model = dryResults$chosen_model,
                         time = j) %>% rbind(compare, .)
   
 }
    
model.result <- data.frame()
# G test to get p value and chosen_model
for(y in unique(compare$gene)){
  
  observed <- as.data.frame(table(compare$chosen_model[compare$gene == y]))
  
  if(nrow(observed) == 1){P = 0}else{
    P <- GTest(observed$Freq, p = rep(1/nrow(observed), nrow(observed)))$p.value}

  if(P < 0.05){
    
    model <- observed$Var1[which.max(observed$Freq)]
    
    model.result <- data.frame(
      gene = y,
      ID = ID,
      chosen_model = model,
      gtest.p.value = P,
      time = max(observed$Freq)) %>% rbind(model.result, .)
  }
  
}
  
  return(model.result)
})

names(r) <- names(compare)
r <- as.data.frame(do.call(rbind, r))
write.table(r, paste0(outdir, clock_gene, '.03.compare.multiple.times.txt'), quote = F, sep = '\t', row.names = F)
