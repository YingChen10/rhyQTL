

rm(list = ls()) 

library(data.table)
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library("dryR"))
suppressMessages(library("lmtest"))

############################# parameters #########################

dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/"
cat <- "cis_QTL/"

args <- commandArgs(trailingOnly = TRUE)
file <- args[1]


logoutdir <- paste0(dir, cat, 'Log/05_compare_among_4_group/')
if (!file.exists(logoutdir)) {dir.create(logoutdir)}

start <- Sys.time()
write.table(start, paste0(logoutdir, file), sep = '\t', quote = F, row.names = F, col.names = F)

female.indir <- paste0(dir, cat, '11_Gender/01_expr_geno_time/Adipose-Visceral_Omentum.female/')
male.indir <- paste0(dir, cat, '11_Gender/01_expr_geno_time/Adipose-Visceral_Omentum.male/')


female <- readRDS(paste0(female.indir, file, '.rds'))
male <- readRDS(paste0(male.indir, file, '.rds'))
gene.intersect <- intersect(names(male), names(female))
all <- lapply(gene.intersect, function(gene.tmp){
  
  
  female.tmp <- female[[gene.tmp]]
  male.tmp <- male[[gene.tmp]]
  id.intersect <- intersect(names(female.tmp), names(male.tmp))
  
  
  gene.result <- lapply(id.intersect, function(id.tmp){
    message(gene.tmp)
    message(id.tmp)
    a <- female.tmp[[id.tmp]] %>% mutate(gender = 'female')
    b <- male.tmp[[id.tmp]] %>% mutate(gender = 'male')
    
    num.a <- as.data.frame(table(a$genotype))
    num.b <- as.data.frame(table(b$genotype))
    num.a.order <- num.a[order(-num.a$Freq),]
    num.b.order <- num.b[order(-num.b$Freq),]
    
    if(paste0(sort(num.a.order$Var1[c(1,2)])[1], '_', sort(num.a.order$Var1[c(1,2)])[2]) == paste0(sort(num.a.order$Var1[c(1,2)])[1], '_', sort(num.a.order$Var1[c(1,2)])[2])){
      
       tmp <- rbind(filter(a, genotype == num.a.order$Var1[1] | genotype == num.a.order$Var1[2]),
                   filter(b, genotype == num.a.order$Var1[1] | genotype == num.a.order$Var1[2]))
      tmp$info <- paste0(tmp$gender, "_", tmp$genotype)
      sample <- as.data.frame(table(tmp$info))
      sample.size <- min(sample$Freq)
      
      random <- data.frame()
      for(info.tmp in unique(sample$Var1)){
        
        now <- filter(tmp, info == info.tmp)
        random <- now[sample(nrow(now), sample.size), ] %>% rbind(random)
      
        }
      
      result.tmp <- as.numeric(drylm(random$x, random$info, random$time)$parameters[1,])
      
    }else{
        
      result.tmp <-NULL
      
    }
    
    return(result.tmp)
    
  })
  
  names(gene.result) <- id.intersect
  gene.result <- as.data.frame(do.call(rbind, gene.result)) %>% rownames_to_column(., var = 'ID')
  names(gene.result) <- c('ID', 
                           paste0(c('mean', 'a', 'b', 'amp', 'relamp', 'phase'), '_female_geno1'),
                           paste0(c('mean', 'a', 'b', 'amp', 'relamp', 'phase'), '_female_geno2'), 
                           paste0(c('mean', 'a', 'b', 'amp', 'relamp', 'phase'), '_male_geno1'),
                           paste0(c('mean', 'a', 'b', 'amp', 'relamp', 'phase'), '_male_geno2'),
                           'chosen_model', 'chosen_model_BICW', 'chosen_model_mean', 'chosen_model_mean_BICW')
  
  gene.result$gene <- gene.tmp

  return(gene.result)
})


all <- as.data.frame(do.call(rbind, all))
write.table(all, paste0(dir, cat, '11_Gender/05_compare_among_4_group/', file), quote = F, sep = '\t', row.names = F)

a <- data.frame(file = file, s = start, e = Sys.time())
write.table(a, paste0(logoutdir, file, '.complete'), sep = '\t', quote = F, row.names = F, col.names = F)
