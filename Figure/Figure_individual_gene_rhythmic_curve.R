
suppressMessages(library(tidyverse))
suppressMessages(library("lmtest"))
suppressMessages(library(gridExtra))
suppressMessages(library(data.table))
suppressMessages(library(grid))
##################### Function #############################

theme <- theme_bw() +
  theme(axis.line = element_blank(),
        axis.text = element_text(color = "black", size = 8),
        #axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) + 
  theme(legend.position = "none")
# get parameters
get_parameter <- function(x, t, period) {
  x <- unlist(x)
  # regression
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
  
  return(data.frame(p = p.val, amp = amp))
}

# Harmonic regression function
harm_reg_plot <- function(x, t, period, expr, genotype_to_plot, sample.size) {
  x <- unlist(x)
  # regression
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
  
  print(p.val)
  print(amp)
  print(phase)
  
  # generat fit curve data
  t.fit=seq(min(t),max(t),0.1)
  c1.f=a*cos(2*pi*t.fit/period)
  s1.f=b*sin(2*pi*t.fit/period)
  
  expr.fit=mu+c1.f+s1.f
  results <- data.frame(t = c(t, t.fit),
                        expr = c(x, expr.fit),
                        hl = c(rep('obs', time = length(t)), rep('fit', time = length(t.fit))))
  
  if(p.val < 0.0001){
    ggplot() +
      geom_point(data = filter(results, hl == 'obs'), aes(x = t, y = expr), size = 0.1, color = 'grey') +
      geom_line(data = filter(results, hl == 'fit'), aes(x = t, y = expr), color = '#4169E1') +
      theme_bw() + theme(panel.grid = element_blank()) + ylim(floor(min(expr)), ceiling(max(expr))) +
      xlab("Time (h)") + ylab("log2 normalized expression") + 
      scale_x_continuous(limits = c(0, 24), breaks = c(0, 6, 12, 18, 24)) +
      ggtitle(paste0(genotype_to_plot, ' n = ', sample.size)) +
      theme(axis.line = element_blank(),
            axis.text = element_text(color = "black", size = 8),
            #axis.text = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            axis.title = element_text(size = 8),
            plot.title = element_text(size = 8))
  }
  else{
    ggplot() +
      geom_point(data = filter(results, hl == 'obs'), aes(x = t, y = expr), size = 0.1, color = 'grey') +
      geom_line(data = filter(results, hl == 'fit'), aes(x = t, y = mean(expr)),color = '#4169E1') +
      theme_bw() + theme(panel.grid = element_blank()) + ylim(floor(min(expr)), ceiling(max(expr))) +
      #ylim(min(x) - 0.1, max(x) + 0.1) +
      xlab("Time (h)") + ylab("log2 normalized expression")+ 
      scale_x_continuous(limits = c(0, 24), breaks = c(0, 6, 12, 18, 24))+
      ggtitle(paste0(genotype_to_plot, ' n = ', sample.size)) +
      theme(axis.line = element_blank(),
            axis.text = element_text(color = "black", size = 8),
            #axis.text = element_blank(),
            panel.grid = element_blank(),
            panel.border = element_rect(color = "black", fill = NA),
            axis.title = element_text(size = 8),
            plot.title = element_text(size = 8))
  }
}

# plot
Plot <- function(tissue, gene, rQTL, period){
  
  # expression profile
  expression <- fread(paste0(dir,'GTEx_nor_expression/',tissue, '.txt'))
  expression <- as.data.frame(expression)
  expression <- expression[,1:(ncol(expression) - 34)]
  expression$EnsemblID <- sapply(strsplit(expression$EnsemblID, '_'), "[", 1)
  names(expression)[2:ncol(expression)] <- sapply(strsplit(names(expression)[2:ncol(expression)], '\\.'), "[", 2)
  expression <- column_to_rownames(expression, var = 'EnsemblID')
  # expression <- na.omit(expression)
  
  # predicted time information
  time <- fread(paste0(dir, 'GTEx_donor_time_science.txt'))
  
  # genotype information
  region <- paste0(unlist(lapply(strsplit(rQTL, "_"), function(x) x[1])), ':', 
                   unlist(lapply(strsplit(rQTL, "_"), function(x) x[2])), '-', 
                   unlist(lapply(strsplit(rQTL, "_"), function(x) x[2])))
  
  
  head <- read.table(paste0(dir, 'head_sub.txt'))
  output <- system(paste0("/workspace/rsrch1/ychen/miniconda3/bin/tabix ", dir, "GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hwe_MAF059.vcf.gz ", region), intern=TRUE)
  snv <- read.table(textConnection(output), header=F, sep="\t")
  names(snv) <- as.character(head)
  snv <- snv[grepl("A|G|C|T", snv$REF, perl = T), ]
  snv <- snv[grepl("A|G|C|T", snv$ALT, perl = T), ]
  snv <- snv[,10:ncol(snv)]
  
  geno <- as.data.frame(t(snv))
  names(geno) <- 'V1'
  geno <- subset(geno, geno$V1 != '.')
  geno$genotype <- as.numeric(substring(geno$V1,1,1)) + as.numeric(substring(geno$V1,3,3))
  geno <- geno[intersect(row.names(geno), names(expression)),]
  
  num <- as.data.frame(table(geno$genotype))

  
  x <- expression[gene, ]
  t <- time$hour[match(names(x), time$SUBJ.ID)]
  p <- harm_reg_plot(x, t, period, x, 'all', length(x))
  
  
  # compare average expression level
  geno$expr <- as.numeric(select(x, row.names(geno)))
  pp <- ggplot(data = geno, aes(x = as.character(genotype), y = expr)) + geom_boxplot(outlier.shape = NA, outlier.colour = NA, color = '#4169E1') + 
    theme_bw() + theme(panel.grid = element_blank()) + ylim(floor(min(x)), ceiling(max(x))) +
    xlab("Genotype") + ylab("log2 normalized expression") + ggtitle(paste0('all n = ', length(geno$expr))) +
    theme(axis.line = element_blank(),
          axis.text = element_text(color = "black", size = 8),
          panel.grid = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          axis.title = element_text(size = 8),
          plot.title = element_text(size = 8))
    
  
  x0 <- select(expression[gene, ], row.names(geno[geno$genotype == 0,]))
  t0 <- time$hour[match(names(x0), time$SUBJ.ID)]
  if(length(x0) > 1){p0 <- harm_reg_plot(x0, t0, period, x, 0, length(x0))}
  
  x1 <- select(expression[gene, ], row.names(geno[geno$genotype == 1,]))
  t1 <- time$hour[match(names(x1), time$SUBJ.ID)]
  if(length(x1) > 1){p1 <- harm_reg_plot(x1, t1, period, x, 1, length(x1))}
  
  x2 <- select(expression[gene, ], row.names(geno[geno$genotype == 2,]))
  t2 <- time$hour[match(names(x2), time$SUBJ.ID)]
  if(length(x2) > 1){p2 <- harm_reg_plot(x2, t2, period, x, 2, length(x2))}
  
  if(nrow(num) == 3){
    pdf(paste0(outdir, gene, '-', rQTL, '-', tissue, '.pdf'), height = 1.6, width = 7.5)
    #grid.arrange(p0, p1, p2, p, pp, ncol = 5)
    print(p0, vp=viewport(.16, .85, x = .15, y = .5))
    print(p1, vp=viewport(.16, .85, x = .3, y = .5))
    print(p2, vp=viewport(.16, .85, x = .45, y = .5))
    print(p, vp=viewport(.16, .85, x = .6, y = .5))
    print(pp, vp=viewport(.16, .85, x = .75, y = .5))
    grid.text(paste0(tissue, ' ', gene, ' ', rQTL), x = 0.5, y = 0.98, just = c("center"), gp = gpar(fontsize = 8))
    }else{
    pdf(paste0(outdir, gene, '-', rQTL, '-', tissue, '.pdf'), height = 1.6, width = 4.5)
    grid.arrange(p0, p1, p2, p, pp, ncol = 4)
    grid.text(paste0(tissue, ' ', gene, ' ', rQTL), x = 0.5, y = 0.98, just = c("center"), gp = gpar(fontsize = 8))
    }
  dev.off()
  
  a <- rbind(
  get_parameter(x, t, period),
  get_parameter(x0, t0, period),
  get_parameter(x1, t1, period),
  get_parameter(x2, t2, period))
  
  return(a)
  
  }
##################################################

################## parameters ###########################
dir <- "/workspace/rsrch1/ychen/Projects/rQTL/"
cat <- "cis_QTL/"



args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
gene <- args[2]
rQTL <- args[3]
out <- args[4]
period <- 24

outdir <- paste0(dir, cat, 'Figures/Individual_gene_rQTL/', out, '/')

if (!file.exists(outdir)) {dir.create(outdir)}

# # 
# tissue <- 'Adipose-Visceral_Omentum'
# gene <- 'ENSG00000130592.15'
# rQTL <- 'chr11_1598677_G_A_b38'

Plot(tissue, gene, rQTL, period)








# period <- 24
# library(grid)
# pdf('./Figures/04_validation/example.pdf', width = 9, height = 1.8)
# # samples with genotype == 0
# x <- select(expression['ENSG00000203747.10', ], row.names(geno[geno$genotype == 0,]))
# t <- time$hour[match(names(x), time$SUBJ.ID)]
# print(harm_reg_plot(x, t, period), vp=viewport(.2, 1, x = .1, y = .5))
# 
# # samples with genotype == 1
# x <- select(expression['ENSG00000203747.10', ], row.names(geno[geno$genotype == 1,]))
# t <- time$hour[match(names(x), time$SUBJ.ID)]
# print(harm_reg_plot(x, t, period), vp=viewport(.2, 1, x = .3, y = .5))
# 
# # samples with genotype == 1
# x <- select(expression['ENSG00000203747.10', ], row.names(geno[geno$genotype == 2,]))
# t <- time$hour[match(names(x), time$SUBJ.ID)]
# print(harm_reg_plot(x, t, period), vp=viewport(.2, 1, x = .5, y = .5))
# 
# # all samples
# x <- select(expression['ENSG00000203747.10', ], row.names(geno))
# t <- time$hour[match(names(x), time$SUBJ.ID)]
# print(harm_reg_plot(x, t, period), vp=viewport(.2, 1, x = .7, y = .5))
# dev.off()
# #results <- rbind(r1, r2, r3) %>% mutate(type = c(rep('all', nrow(r1)), rep('0', nrow(r2)), rep('1', nrow(r3))))
# 
# 
# ###########################
# # NAc tissue
# ###########################
# dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/human_brain_circadian_RNA/rQTL/"
# cat <- "rQTL_GWASCatalog/"
# tissue <- 'NAc'
# meta <- read.table(paste0(dir, 'meta_mod.txt'), header = T)
# 
# # expression profile
# expression <- fread(paste0(dir, 'GTEx_nor_expression/', tissue, '.txt'))
# expression <- na.omit(expression)
# expression <- data.frame(expression, check.names = FALSE)
# expression <- column_to_rownames(expression, var = 'V1')
# 
# 
# # input time
# time <- fread(paste0(dir, 'meta_data_github.txt'))
# time$hour <- time$CorrectedTOD + 6
# names(time)[names(time) == 'pair'] <- 'SUBJ.ID'
# head(time)
# 
# # genotype information
# snv <- read.table('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/human_brain_circadian_RNA/vcf_combine/filter/NAc_GWAScat.vcf')
# head <- read.table('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/human_brain_circadian_RNA/vcf_combine/filter/NAc_GWAScat.head')
# names(snv) <- c('CHR', as.character(head$V1))
# snv <- snv[grepl("A|G|C|T", snv$REF, perl = T), ]
# snv <- snv[grepl("A|G|C|T", snv$ALT, perl = T), ]
# 
# meta[match(names(snv)[10:ncol(snv)], meta$Sample_geo_accession),]$Sample -> names(snv)[10:ncol(snv)]
# snv <- select(snv, ID, 10:ncol(snv))
# genotype <- as.data.frame(t(apply(snv, 1, function(x){str_split_fixed(x, ':', 8)[ ,1]})))
# names(genotype) <- names(snv)
# region
# geno <- as.data.frame(t(genotype[genotype$ID == 'chr2_61337114_A_C_b38', -1]))
# names(geno) <- 'V1'
# geno <- subset(geno, geno$V1 != './.')
# geno$genotype <- as.numeric(substring(geno$V1,1,1)) + as.numeric(substring(geno$V1,3,3))
# geno <- geno[intersect(row.names(geno), names(expression)),]
# 
# period <- 24
# library(grid)
# pdf('./Figures/04_validation/NAc_example.pdf', width = 9, height = 1.8)
# # samples with genotype == 0
# x <- select(expression['ENSG00000203747', ], row.names(geno[geno$genotype == 0,]))
# t <- time$hour[match(names(x), time$SUBJ.ID)]
# print(harm_reg_plot(x, t, period), vp=viewport(.2, 1, x = .1, y = .5))
# 
# # samples with genotype == 1
# x <- select(expression['ENSG00000203747', ], row.names(geno[geno$genotype == 1,]))
# t <- time$hour[match(names(x), time$SUBJ.ID)]
# print(harm_reg_plot(x, t, period), vp=viewport(.2, 1, x = .3, y = .5))
# 
# # samples with genotype == 1
# x <- select(expression['ENSG00000203747', ], row.names(geno[geno$genotype == 2,]))
# t <- time$hour[match(names(x), time$SUBJ.ID)]
# print(harm_reg_plot(x, t, period), vp=viewport(.2, 1, x = .5, y = .5))
# 
# # all samples
# x <- select(expression['ENSG00000203747', ], row.names(geno))
# t <- time$hour[match(names(x), time$SUBJ.ID)]
# print(harm_reg_plot(x, t, period), vp=viewport(.2, 1, x = .7, y = .5))
# dev.off()
# 
# # Harmonic regression function
# harm_reg_plot <- function(x, t, period){
#   x <- unlist(x)
#   # regression
#   n=length(x)
#   fit0=lm(x~1)
#   c=cos(2*pi*t/period)
#   s=sin(2*pi*t/period)
#   fit1=lm(x~c+s)
#   mu=coef(fit1)[1]
#   a=coef(fit1)[2]
#   b=coef(fit1)[3]
#   p.val=lrtest(fit1, fit0)$Pr[2]
#   amp=2*sqrt(a^2+b^2)
#   phase=atan2(b,a)%%(2*pi)
#   phase=period*phase/(2*pi)
#   
#   # generate fit curve data
#   t.fit=seq(min(t),max(t),0.1)
#   c1.f=a*cos(2*pi*t.fit/period)
#   s1.f=b*sin(2*pi*t.fit/period)
#   
#   expr.fit=mu+c1.f+s1.f
#   results <- data.frame(t = c(t, t.fit),
#                         expr = c(x, expr.fit),
#                         hl = c(rep('obs', time = length(t)), rep('fit', time = length(t.fit))))
#   #return(results)
#   
#   if(p.val < 0.06){
#     ggplot() +
#       geom_point(data = filter(results, hl == 'obs'), aes(x = t, y = expr), size = 0.2) +
#       geom_line(data = filter(results, hl == 'fit'), aes(x = t, y = expr)) +
#       theme_bw() + theme(panel.grid = element_blank()) + 
#       xlab("Time (h)") + ylab("log2 normalized expression")+
#       scale_y_continuous(limits = c(0, 4.5), breaks = c(0, 1, 2, 3, 4)) +
#       scale_x_continuous(limits = c(0, 24), breaks = c(0, 6, 12, 18, 24))
#   }
#   else{
#     ggplot() +
#       geom_point(data = filter(results, hl == 'obs'), aes(x = t, y = expr), size = 0.2) +
#       geom_line(data = filter(results, hl == 'fit'), aes(x = t, y = mean(expr))) +
#       theme_bw() + theme(panel.grid = element_blank()) + 
#       xlab("Time (h)") + ylab("log2 normalized expression")+
#       scale_y_continuous(limits = c(0, 4.5), breaks = c(0, 1, 2, 3, 4)) + 
#       scale_x_continuous(limits = c(0, 24), breaks = c(0, 6, 12, 18, 24))
#   }
# }
# 
# library(tidyverse)
# x <- as.data.frame(t(expression['ENSG00000203747.10', ]))
# x <- rownames_to_column(x, var = 'id')
# names(x) <- c('id', 'expr')
# 
# 
# geno <- as.data.frame(t(snv))
# names(geno) <- 'V1'
# geno <- subset(geno, geno$V1 != '.')
# geno$genotype <- as.numeric(substring(geno$V1,1,1)) + as.numeric(substring(geno$V1,3,3))
# geno <- geno[intersect(row.names(geno), names(expression)),]
# geno <- rownames_to_column(geno, var = 'id')
# geno1 <- merge(geno, x, by = 'id')
# 
# geno2 <- geno1
# geno2$genotype <- 3
# geno2 <- rbind(geno1, geno2)
# head(geno1)
# 
# p <- ggplot(data = geno2, aes(x = as.character(genotype), y = expr)) + geom_boxplot(outlier.shape = NA) + geom_jitter(size = 0.2) +
#   theme_bw() + theme(panel.grid = element_blank())
# p
# library(grid)
# pdf('./Figures/04_validation/NAc_example_eqtl.pdf',  width = 1.8, height = 1.8)
# print(p, vp=viewport(1, 1, x = .5, y = .5))
# dev.off()
# 
# t.test(geno1$expr[geno1$genotype == 0], geno1$expr[geno1$genotype == 1])
# t.test(geno1$expr[geno1$genotype == 2], geno1$expr[geno1$genotype == 1])
# t.test(geno1$expr[geno1$genotype == 2], geno1$expr[geno1$genotype == 0])
# 
# wilcox.test(geno1$expr[geno1$genotype == 0], geno1$expr[geno1$genotype == 1])
# wilcox.test(geno1$expr[geno1$genotype == 2], geno1$expr[geno1$genotype == 1])
# wilcox.test(geno1$expr[geno1$genotype == 2], geno1$expr[geno1$genotype == 0])

