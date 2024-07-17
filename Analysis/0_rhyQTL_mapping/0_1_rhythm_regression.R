
suppressMessages(library(data.table))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library("dryR"))
suppressMessages(library("lmtest"))

############################# parameters #########################
dir <- "Projects/rQTL/"
cat <- "cis_QTL/00_rQTL_mapping/00_Genotype/"
outdir <- paste0(dir, 'cis_QTL/00_rQTL_mapping/01_Rhythm_regression_filter/')
logoutdir <-  paste0(dir, 'cis_QTL/Log/01_Rhythm_regression_filter/')

args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
pos_file <- args[2]
# 
# tissue <- 'Vagina'
# pos_file <- 'HLA_split_pos_ac'


if (!file.exists(paste0(outdir, tissue))) {
  dir.create(paste0(outdir, tissue))}

if (!file.exists(paste0(logoutdir, tissue))) {
  dir.create(paste0(logoutdir, tissue))}
a <- Sys.time()
write.table(a, paste0(logoutdir, tissue, '/', pos_file), quote = F, row.names = F, col.names = F)

###################### Function ###################################
# Harmonic regression function
harm_reg <- function(x, t, period){
  n=length(x)
  fit0=lm(x~1)
  c=cos(2*pi*t/period)
  s=sin(2*pi*t/period)
  fit1=lm(x~c+s)
  mu=coef(fit1)[1]
  a=coef(fit1)[2]
  b=coef(fit1)[3]
  p.val=lrtest(fit1, fit0)$Pr[2]
  amp=2*sqrt(a^2+b^2)
  phase=atan2(b,a)%%(2*pi)
  phase=period*phase/(2*pi)
  return(c(pval=p.val,phase=phase,amp=amp,mu=mu,a=a,b=b, period=period, R2=summary(fit1)$r.squared))
}


# Harmonic regression for each genotype group
Regression <- function(all_results, time, expression){
  
  geno_regression <- list()
  
  for(i in 1:length(all_results)){
    
    gene_list <- all_results[[i]]
    gene <- names(all_results)[i]
    
    regression.results <- lapply(gene_list, function(genotype_list){
      
      geno_tmp = genotype_list$geno
      tmp <- data.frame(genotype = geno_tmp$genotype,
                        time = time$hour[match(rownames(geno_tmp), time$SUBJ.ID)],
                        x = as.numeric(select(filter(expression, EnsemblID == gene), row.names(geno_tmp))))
      
      # calculate sample size for each genotype subgroup and sort genotype subgroup anccording to sample (from large to small) 
      num <- as.data.frame(table(factor(tmp$genotype, levels = 0:2)))
      num <- num[order(-num$Freq), ]
      
      # fit rhythmicity in each genotype subgroup
      regression <- num$Freq
      for(j in 1:2){
        
        tmp_geno <- filter(tmp, genotype == num$Var1[j])
        as.numeric(harm_reg(x = tmp_geno$x, t = tmp_geno$time, period = period))
        regression <- c(regression, as.numeric(harm_reg(x = tmp_geno$x, t = tmp_geno$time, period = period)))
        
      }
      
      return(regression)
      
    })
    
    regression <- as.data.frame(do.call(rbind,  regression.results))
    names(regression) <- c('smaple.size.0', 'smaple.size.1', 'smaple.size.2', paste0(name, '_0'), paste0(name, '_1'))

    regression$qval_0=p.adjust(regression$pval_0, "BH")
    regression$qval_1=p.adjust(regression$pval_1, "BH")
    
    filter_id <- intersect(row.names(regression[rowSums(select(regression, contains('pval')) < 0.01) > 0,]), row.names(regression[rowSums(select(regression, contains('amp.c')) > log2(1.5)) > 0,]))
    regression <- regression[filter_id,]
    
if(nrow(regression) != 0){geno_regression[[gene]] <- regression}
    
  }
  
  return(geno_regression)
}
######################################################################


###### import files ##############
### import time information
time <- fread(paste0(dir, 'GTEx_donor_time_science.txt'))

# import expression file
expression <- fread(paste0(dir,'GTEx_nor_expression/',tissue, '.txt'))
expression <- as.data.frame(expression)
expression <- expression[,1:(ncol(expression) - 34)]
expression$EnsemblID <- sapply(strsplit(expression$EnsemblID, '_'), "[", 1)
names(expression)[2:ncol(expression)] <- sapply(strsplit(names(expression)[2:ncol(expression)], '\\.'), "[", 2)

# import genotype information
geno <- readRDS(paste0(dir, cat, tissue, '/', pos_file, '.rds'))

period <- 24
name <- c("pval","phase.s","amp.c","mu.(Intercept)","a.c","b.s","period","R2")
########################################################## Main ###############################################

######################### 1. Regression ################################
filter.result <- Regression(geno, time, expression)
saveRDS(filter.result, paste0(outdir, tissue, '/', pos_file, '.rds'))
###############


filter.geno.list <- c()
for(i in 1:length(filter.result)){

  
  gene <- names(filter.result)[i]
  fit.filter <- filter.result[[gene]]
  
  if(nrow(fit.filter) == 0){next}
  
  filter.geno <- lapply(row.names(fit.filter), function(id){
    
    dif.rhy <- c()
    
    # get genotype information for current SNP in the gene
    tmp_geno <- geno[[gene]][[id]]$geno
    
    return(tmp_geno)
    })
  names(filter.geno) <- row.names(fit.filter)
  
  filter.geno.list[[gene]] <- filter.geno
}

saveRDS(filter.geno.list, paste0(outdir, tissue, '/', pos_file, '.geno.rds'))

b <- data.frame(tissue = tissue, file = pos_file, s = a, e = Sys.time(), gene.filter = length(filter.geno.list))
write.table(b, paste0(logoutdir, tissue, '/', pos_file, '.final'), quote = F, row.names = F, col.names = F)

