
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


## filter SNV and generate genotype list and down sampling
Filter_SNV <- function(tmp_snv, min_size){
  
  # calculate sample size in each genotype group
  geno <- as.data.frame(t(tmp_snv))
  names(geno) <- 'V1'
  geno <- subset(geno, geno$V1 != '.')
  geno$genotype <- as.numeric(substring(geno$V1,1,1)) + as.numeric(substring(geno$V1,3,3))
  num <- as.data.frame(table(geno$genotype))
  
  # only keep wt group with > min size and keep 1 group of mut group (if only 1 mutation group > 50 keep this group if 2 mutation groups > 50, keep the group with larger sample size)
  num_cp <- filter(num, Freq >= min_size)
  
  if(sum(num_cp$Var1 == 0) > 0 & nrow(num_cp) > 1){
    num_mut <- filter(num_cp, Var1 != 0)
    if(nrow(num_mut) == 1){mut_group = num_mut$Var1}else{
      # only keep mutation group with larger sample size
      if(num_mut$Freq[1] == num_mut$Freq[2]){mut_group = 1}else{
        mut_group = num_mut$Var1[num_mut$Freq == max(num_mut$Freq)]}}
    
    # subset geno and do down sampling
    geno <- filter(geno, genotype == 0 | genotype == mut_group)
    geno <- DownSample(geno)
    geno$genotype[geno$genotype != 0] <- 1
    
    # return
    result <- list(geno = geno, num = num, size2test = nrow(geno)/2)
    
    return(result)
  }
  else{
    
    return(NULL)
  }
}

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

source('~/Projects/Project03_human_circadian/Scripts/rQTL/cis_rQTL/Script/05_trans_rGene/dff_rhythm.R')

############################## Main #######################################
### input time information
time <- fread(paste0(dir, 'GTEx_donor_time_science.txt'))

### head of snv file
head <- read.table(paste0(dir, 'head_sub.txt'))

###### input snv file
# get snv id
rqtl <- read.table(paste0(dir, cat, '/13_Results/00_all_rhyQTL/', tissue, '.txt'), header = T)
snv_id <- filter(rqtl, name == clock_gene)


# get all indivadual's snv information
chr <- unique(str_split_fixed(snv_id$ID, '_', 5)[ ,1])
pos <- as.numeric(str_split_fixed(snv_id$ID, '_', 5)[ ,2])
start <- pos[order(pos)][1]
end <- pos[order(pos)][length(pos)]
region <- paste0(chr, ':', start, '-', end)

# retrieve SNVs in TSS-1M flanking region  for each gene
output <- system(paste0("/workspace/rsrch1/ychen/miniconda3/bin/tabix ", dir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hwe_MAF059.vcf.gz ", region), intern=TRUE)
snv <- read.table(textConnection(output), header=F, sep="\t")
names(snv) <- as.character(head)
snv <- snv[grepl("A|G|C|T", snv$REF, perl = T), ]
snv <- snv[grepl("A|G|C|T", snv$ALT, perl = T), ]
snv <- filter(snv, ID %in% snv_id$ID)


# input expression file
expression <- fread(paste0(dir,'GTEx_nor_expression/', tissue, '.txt'))
expression <- as.data.frame(expression)
expression <- expression[,1:(ncol(expression) - 34)]
expression$EnsemblID <- sapply(strsplit(expression$EnsemblID, '_'), "[", 1)
names(expression)[2:ncol(expression)] <- sapply(strsplit(names(expression)[2:ncol(expression)], '\\.'), "[", 2)
expression <- column_to_rownames(expression, var = 'EnsemblID')


print('rQTL mapping statring!')

r <- lapply(snv$ID, function(x){
  
  # print(x)
  
  # get information for the SNP ID and subset samples with both expression and genotype infomation
  tmp_snv <- filter(snv, ID==x) %>% select(., intersect(names(expression), names(snv)))
  tmp_snv <- Filter_SNV(tmp_snv, min_size)
  
  if (is.null(tmp_snv)){return(NULL)}
  geno <- tmp_snv$geno
  
  tmp <- c()
  # Harmonic regression for each genotype
  for(i in 0:1){
    
    geno_tmp <- filter(geno, genotype == i)
    expression_tmp <- select(expression, row.names(geno_tmp))
    time.point <- time$hour[match(row.names(geno_tmp), time$SUBJ.ID)]
    
    assign(paste0('data.fit.', i), lapply(rownames(expression_tmp), function(x) harm_reg(as.numeric(expression_tmp[x,]), time.point, period)))
    
    assign(paste0('t', i), time.point %>% unlist())
    assign(paste0('e', i), expression_tmp %>% as.matrix())
    
  }
  
  HANOVA.no.norm  = as.numeric(HANOVA(t(e0), t(e1), t0, t1, 24, norm = FALSE)$p.value)
  HANOVA.norm  = as.numeric(HANOVA(t(e0), t(e1), t0, t1, 24, norm = TRUE)$p.value)
  
  #### integrate results for all genotyoes
  names(data.fit.0) <- row.names(expression)
  data.fit.0 <- as.data.frame(do.call(rbind, data.fit.0))
  data.fit.0$qval <- p.adjust(data.fit.0$pval, "BH")
  
  names(data.fit.1) <- row.names(expression)
  data.fit.1 <- as.data.frame(do.call(rbind, data.fit.1))
  data.fit.1$qval <- p.adjust(data.fit.1$pval, "BH")
  
  results <- cbind(data.fit.0, data.fit.1)
  message(nrow(results))
  names(results) <- c(paste0(names(data.fit.0), '_0'), paste0(names(data.fit.0), '_1'))
  results$max_amp <- apply(select(results, contains("amp.c")), 1, function(x) max(x))
  results$max_amp_idx <- apply(select(results, contains("amp.c_")), 1, function(x) which.max(x))
  results$pval <- apply(select(results, contains("pval_"), max_amp_idx), 1, function(x) (x[x[3]]))
  results$qval <- apply(select(results, contains("qval_"), max_amp_idx), 1, function(x) (x[x[3]]))
  results$HANOVA.no.norm <- HANOVA.no.norm
  results$HANOVA.norm <- HANOVA.norm
  results$HANOVA.norm.qval <- p.adjust(results$HANOVA.norm, "BH")
  results <- filter(results, pval <= 0.01 & max_amp > 0.5)
  
  # ### get significant results and Run dryR
  # dryR_test <- intersect(row.names(results[rowSums(select(results, contains('pval')) < 0.01) > 0,]), row.names(results[rowSums(select(results, contains('amp.c')) > 0.5) > 0,]))
  # 
  # if (length(dryR_test) > 0){
  #   
  #   expression.test <- na.omit(expression[dryR_test,])
  #   dryList = drylm(select(expression.test, row.names(geno)), as.character(geno$genotype), time$hour[match(row.names(geno), time$SUBJ.ID)])
  #   
  #   dryResults <- dryList[["parameters"]] # coefficients: phase, amplitude and mean for each group
  #   dryResults <- filter(dryResults, chosen_model != 1 & chosen_model != 4)
  #   
  #   subset <- results[row.names(dryResults), ]
  #   dryResults <- cbind(subset, dryResults)
  #   dryResults$sample.size_test = tmp_snv$size2test
  # }else{
  #   
  #   dryResults <- NULL
  # }
  # 
  return.results <- list(tmp_snv = tmp_snv, dryResults = results)
  return(return.results)
})

names(r) <- snv$ID

sample <- lapply(r, function(x){
  return(x$tmp_snv)
})

sig <- lapply(r, function(x){
  return(x$dryResults)
})

print('output results!')

saveRDS(sample, paste0(outdir, clock_gene, '.01.sample.rds'))
saveRDS(sig, paste0(outdir, clock_gene, '.01.sig.rds'))


b <- as.data.frame(do.call(rbind, sig))
b <- rownames_to_column(b, var = 'ID')
b$gene <- str_split_fixed(b$ID, 'b38.', 2)[ ,2]
b$ID <- paste0(str_split_fixed(b$ID, 'b38.', 2)[ ,1], 'b38')
annot <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/gencode.v26.GRCh38.genes.annot')
annot <- filter(annot, V2 == 'gene')[,6:7]
names(annot) <- c('gene', 'name')
all <- merge(annot, b, by = 'gene')
write.table(all, paste0(outdir, clock_gene, '.01.regression.txt'), row.names = F, sep = '\t', quote = F)


