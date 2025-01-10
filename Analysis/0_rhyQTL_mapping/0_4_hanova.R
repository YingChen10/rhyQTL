
rm(list=ls())

suppressWarnings({
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
})


HANOVA <- function (val1, val2,
                    times1, times2,
                    period, norm = TRUE){
  
  #prepare val1 and val2 if they contain only one time series
  if(length(dim(val1)) == 1){
    val1 <- matrix(val1, ncol=1)
  }
  if(length(dim(val2)) == 1){
    val2 <- matrix(val2, ncol=1)
  }
  
  #normalisation
  if(norm){
    val1 <- apply(val1, 2, function(series){
      normval <- series/mean(series, na.rm=TRUE)
      return(normval)
    })
    
    val2 <- apply(val2, 2, function(series){
      normval <- series/mean(series, na.rm=TRUE)
      return(normval)
    })
  }
  
  #prepare matrices for Rfit calculations
  timevec <- c(times1, times2)/period*2*pi
  
  #vectors for combined fit
  cosB <- cos(timevec)
  sinB <- sin(timevec)
  
  #vetors for individual Fits
  cosD <- cosB*c(rep(0,length(times1)),rep(1,length(times2)))
  sinD <- sinB*c(rep(0,length(times1)),rep(1,length(times2)))
  
  
  #combine both value matrices
  valC <- rbind(val1, val2)
  
  tp1 <- length(times1)
  tp2 <- length(times2)
  
  muchna <- apply(val1, 2, function(x) sum(is.na(x))-tp1 >= -3) |
    apply(val2, 2, function(x)  sum(is.na(x))-tp1 >= -3)
  
  xred <- cbind(cosB, sinB)
  xfull <- cbind(cosB, sinB, cosD, sinD)
  
  relist <- function(obj){
    len <- ncol(obj$coefficients)
    lm.list <- lapply(seq_len(len), function(index){
      ret <- list(coefficients = obj$coefficients[,index],
                  residuals = obj$residuals[,index],
                  effects = obj$effects[,index],
                  df.residual = obj$df.residual,
                  model = obj$model,
                  call = obj$call)
      class(ret) = 'lm'
      ret <- as(ret, 'lm')
      return(ret)
    })
    
  }
  
  
  #apply fitting
  has.na <- TRUE#!all(!is.na(valC))
  
  
  results <- do.call(rbind,
                     mclapply(seq_len(ncol(val1)), function(index){
                       vals <- valC[, index]
                       if(muchna[index]){
                         retlist <- c(
                           NA,
                           NA,
                           NA)
                         return(retlist)
                       }
                       f.r <- lm(vals ~ xred, na.action = na.omit)
                       f.f <- lm(vals ~ xfull, na.action = na.omit)
                       
                       dt1 <- anova(f.r, f.f)
                       
                       if("try-error" %in% class(dt1)){
                         dt1 <- list('Pr(>F)'=1, F=0)
                         diff <- 0
                       } else{
                         diff <- sqrt(sum(coefficients(f.f)[4:5]^2))
                       }
                       retlist <- c(
                         dt1$'Pr(>F)'[2],
                         dt1$F[2],
                         diff)
                       return(retlist)
                     }, mc.preschedule=TRUE, mc.cleanup=TRUE))
  
  resultdf <- data.frame(results)
  colnames(resultdf) <- c('p.value', 'F', 'diff')
  return(resultdf)
  
}

args <- commandArgs(trailingOnly = T)
tissue <- args[1]
pos_file <- args[2]



dir <- '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/'
cat <- 'cis_QTL/'
outdir <- paste0(dir, cat, '00_rQTL_mapping/04_hanova/')
logoutdir <- paste0(dir, 'cis_QTL/Log/04_hanova/')

if (!file.exists(outdir)) {dir.create(outdir, recursive = T)}
if (!file.exists(logoutdir)) {dir.create(logoutdir, recursive = T)}
if (!file.exists(paste0(outdir, tissue))) {dir.create(paste0(outdir, tissue))}
if (!file.exists(paste0(logoutdir, tissue))) {dir.create(paste0(logoutdir, tissue))}

a <- Sys.time()
write.table(a, paste0(logoutdir, tissue, '/', pos_file), quote = F, row.names = F, col.names = F)

###### import files ##############
### import time information
time <- fread(paste0(dir, 'GTEx_donor_time_science.txt'))

# import expression file
expression <- read.table(paste0(dir,'00_data/CPM_covariate_remove/', tissue, '.txt'), header = T)
names(expression) <- str_split_fixed(names(expression), '\\.', 2)[ ,2]
expression$EnsemblID <- str_split_fixed(row.names(expression), '_', 2)[ ,1]

# import genotype information
geno <- readRDS(paste0(dir, cat, '00_rQTL_mapping/00_Genotype/', tissue, '/', pos_file, '.rds'))


hanova <- lapply(1:length(geno), function(i){
  
  gene_list <- geno[[i]]
  gene <- names(geno)[i]
  
  hanova.results <- lapply(gene_list, function(genotype_list){
    
    geno_tmp <- genotype_list$geno
    tmp <- data.frame(genotype = geno_tmp$genotype,
                      time = time$hour[match(rownames(geno_tmp), time$SUBJ.ID)],
                      x = as.numeric(select(filter(expression, EnsemblID == gene), row.names(geno_tmp))))
    
    # calculate sample size for each genotype and sort genotype anccording to sample (from large to samll) 
    num <- as.data.frame(table(factor(tmp$genotype, levels = 0:2)))
    num <- num[order(-num$Freq), ]
    
    # compare and get p value
    t1 = tmp %>% filter(genotype == num$Var1[1]) %>% select(time) %>% unlist()
    e1 = tmp %>% filter(genotype == num$Var1[1]) %>% select(x) %>% as.matrix()
    t2 = tmp %>% filter(genotype == num$Var1[2]) %>% select(time) %>% unlist()
    e2 = tmp %>% filter(genotype == num$Var1[2]) %>% select(x) %>% as.matrix()
   
    tmp.result <- data.frame(
      
      geno.0 = num$Var1[1],
      geno.1 = num$Var1[2],
      geno.2 = num$Var1[3],
      geno.0.num = num$Freq[1],
      geno.1.num = num$Freq[2],
      geno.2.num = num$Freq[3],
      HANOVA.norm  = as.numeric(HANOVA(e1, e2, t1, t2, 24, norm = TRUE)$p.value))
    
    return(tmp.result)
    
  })
  
  names(hanova.results) <- names(gene_list)
  hanova.results <- as.data.frame(do.call(rbind, hanova.results)) %>% rownames_to_column(., var = 'ID')

  return(hanova.results)
})

names(hanova) <- names(geno)
saveRDS(hanova, paste0(outdir, tissue, '/', pos_file, '.rds'))

b <- data.frame(tissue = tissue, file = pos_file, s = a, e = Sys.time())
write.table(b, paste0(logoutdir,  tissue, '/', pos_file, '.complete'), sep = '\t', quote = F, row.names = F, col.names = F)

