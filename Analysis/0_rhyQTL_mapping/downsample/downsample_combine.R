setwd('F:/Projects/Project03_human_circadian/rQTL/cis_QTL/11_downsample/')


tissue.list <- c('Adipose-Visceral_Omentum_subsample.150',
                 'Adipose-Visceral_Omentum_subsample.200',
                 'Adipose-Visceral_Omentum_subsample.250',
                 'Adipose-Visceral_Omentum_subsample.300',
                 'Adipose-Visceral_Omentum_subsample.350',
                 'Adipose-Visceral_Omentum_subsample.400',
                 'Adipose-Visceral_Omentum_subsample.450')

file.list <- list.files('./12_hanova/Adipose-Visceral_Omentum_subsample.150/') %>% gsub('.rds', '', .)

all <- lapply(tissue.list, function(tissue){
  
  
  file.tmp <- lapply(file.list, function(file){
    
    # regerssion
    regression.tmp <- readRDS(paste0('01_Rhythm_regression/', tissue, '/', file, '.rds'))
    
    regression.tmp <- lapply(regression.tmp, function(tmp){
      
      tmp <- rownames_to_column(tmp, var = 'ID')
      return(tmp)
      
    })
    
    regression.tmp <- as.data.frame(do.call(rbind,  regression.tmp)) %>% rownames_to_column(., var = 'gene')
    regression.tmp$gene <- paste0(str_split_fixed(regression.tmp$gene, '\\.', 3)[ ,1], '.', str_split_fixed(regression.tmp$gene, '\\.', 3)[ ,2])
    regression.tmp$max_amp <- apply(select(regression.tmp, contains("amp.c")), 1, function(x) max(x))
    regression.tmp$max_amp_idx <- apply(select(regression.tmp, contains("amp.c_")), 1, function(x) which.max(x))
    regression.tmp$pval <- apply(select(regression.tmp, contains("pval_"), max_amp_idx), 1, function(x) (x[x[3]]))
    regression.tmp$qval <- apply(select(regression.tmp, contains("qval_"), max_amp_idx), 1, function(x) (x[x[3]]))
    regression.tmp <- filter(regression.tmp, pval < 5e-4)
    
    # hanaove
    hanova.tmp <- readRDS(paste0('04_hanova/', tissue, '/', file, '.rds'))
    hanova.tmp <- as.data.frame(do.call(rbind,  hanova.tmp))
    hanova.tmp <- hanova.tmp[,c(1:4, 8)]
    hanova.tmp <- rownames_to_column(hanova.tmp, var = 'gene')
    hanova.tmp$gene <- paste0(str_split_fixed(hanova.tmp$gene, '\\.', 3)[ ,1], '.', str_split_fixed(hanova.tmp$gene, '\\.', 3)[ ,2])
    
    # combine
    combine <- merge(regression.tmp, hanova.tmp, by = c('gene', 'ID'))
    
    return(combine)
    
  })
    
  file.tmp <- as.data.frame(do.call(rbind, file.tmp)) 
  file.tmp$downsample.szie <- str_split_fixed(tissue, '\\.', 2)[ ,2]
  
  return(file.tmp)
    
  })
names(all) <- tissue.list
saveRDS(all, 'regreseion_hanova_combine.rds')  


############################################################
#
# compare all sample
#
###########################################################
compare.all <- lapply(tissue.list, function(tissue){
  
  file.tmp <- lapply(file.list, function(file){
  
  # compare all samples 
  compare.tmp <- fread(paste0('02_Rhythm_compare_all/', tissue, '/', file))
  
  return(compare.tmp)
  
  })
  
  file.tmp <- as.data.frame(do.call(rbind, file.tmp)) 
  file.tmp$downsample.szie <- str_split_fixed(tissue, '\\.', 2)[ ,2]
  
  return(file.tmp)

})
names(compare.all) <- tissue.list  

count <- lapply(1:length(compare.all), function(i){
  
  compare.tmp <- compare.all[[i]]
  regression.tmp <- all[[i]]
  names(compare.tmp)[names(compare.tmp) == 'Gene'] <- 'gene'
  
  tmp <- merge(regression.tmp, compare.tmp, by = c('gene', 'ID'))
  rhyGene.num <- length(unique(tmp$gene))
  rhyQTL.num <- length(unique(tmp$ID))
  tmp$p_adj <- p.adjust(tmp$HANOVA.norm, method = "BH")
  tmp <- filter(tmp, p_adj < 0.01)
  
  return(
    
    data.frame(random.sample.size = unique(tmp$downsample.szie.x),
               rhyGene.num = rhyGene.num,
               rhyQTL.num = rhyQTL.num)
    
  )
  
})

count <- as.data.frame(do.call(rbind, count))
write.table(count, 'compare_all_sample_count.txt', quote = F, sep = '\t', row.names = F)

############################################################
#
# compare downsample one time
#
###########################################################
compare.all <- lapply(tissue.list, function(tissue){
  
  file.tmp <- lapply(file.list, function(file){
    
    # compare all samples 
    compare.tmp <- fread(paste0('02_Rhythm_compare/', tissue, '/', file))
    
    return(compare.tmp)
    
  })
  
  file.tmp <- as.data.frame(do.call(rbind, file.tmp)) 
  file.tmp$downsample.szie <- str_split_fixed(tissue, '\\.', 2)[ ,2]
  
  return(file.tmp)
  
})
names(compare.all) <- tissue.list  

compare.downsample.count <- lapply(1:length(compare.all), function(i){
  
  compare.tmp <- compare.all[[i]]
  regression.tmp <- all[[i]]
  names(compare.tmp)[names(compare.tmp) == 'Gene'] <- 'gene'
  
  tmp <- merge(regression.tmp, compare.tmp, by = c('gene', 'ID'))
  rhyGene.num <- length(unique(tmp$gene))
  rhyQTL.num <- length(unique(tmp$ID))
  tmp$p_adj <- p.adjust(tmp$HANOVA.norm, method = "BH")
  tmp <- filter(tmp, p_adj < 0.01)
  
  return(
    
    data.frame(random.sample.size = unique(tmp$downsample.szie.x),
               rhyGene.num = rhyGene.num,
               rhyQTL.num = rhyQTL.num)
    
  )
  
})

compare.downsample.count <- as.data.frame(do.call(rbind, compare.downsample.count))
write.table(compare.downsample.count, 'compare_downsample_1time_count.txt', quote = F, sep = '\t', row.names = F)

############################################################
#
# compare downsample 20 times
#
###########################################################
compare.all <- lapply(tissue.list, function(tissue){
  
  file.tmp <- lapply(file.list, function(file){
    
    # compare all samples 
    compare.tmp <- fread(paste0('03_Rhythm_compare_multiple_times/', tissue, '/', file))
    
    return(compare.tmp)
    
  })
  
  file.tmp <- as.data.frame(do.call(rbind, file.tmp)) 
  file.tmp$downsample.szie <- str_split_fixed(tissue, '\\.', 2)[ ,2]
  
  return(file.tmp)
  
})
names(compare.all) <- tissue.list  

compare.downsample.count <- lapply(1:length(compare.all), function(i){
  
  compare.tmp <- compare.all[[i]]
  regression.tmp <- all[[i]]
  names(compare.tmp)[names(compare.tmp) == 'Gene'] <- 'gene'
  
  tmp <- merge(regression.tmp, compare.tmp, by = c('gene', 'ID'))
  tmp <- filter(tmp, gtest.p.value < 0.05)
  rhyGene.num <- length(unique(tmp$gene))
  rhyQTL.num <- length(unique(tmp$ID))
  tmp$p_adj <- p.adjust(tmp$HANOVA.norm, method = "BH")
  tmp <- filter(tmp, p_adj < 0.01)
  
  return(
    
    data.frame(random.sample.size = unique(tmp$downsample.szie.x),
               rhyGene.num = rhyGene.num,
               rhyQTL.num = rhyQTL.num)
    
  )
  
})

compare.downsample.20time.count <- as.data.frame(do.call(rbind, compare.downsample.count))
write.table(compare.downsample.20time.count, 'compare_downsample_20time_count.txt', quote = F, sep = '\t', row.names = F)

#########################
#
# combine and polt
#
#########################
setwd('F:/Projects/Project03_human_circadian/rQTL/cis_QTL/11_downsample/')
a <- fread('compare_all_sample_count.txt') %>% rbind(., data.frame(random.sample.size = 100, rhyGene.num = 0, rhyQTL.num = 0, a = 0, b = 0))
b <- fread('compare_downsample_1time_count.txt')%>% rbind(., data.frame(random.sample.size = 100, rhyGene.num = 0, rhyQTL.num = 0, a = 0, b = 0))
c <- fread('compare_downsample_20time_count.txt')%>% rbind(., data.frame(random.sample.size = 100, rhyGene.num = 0, rhyQTL.num = 0, a = 0, b = 0))

theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) 

p1 <- ggplot(data = a, aes(x = random.sample.size, y = rhyGene.num)) + 
  geom_point() +
  geom_line(color = 'blue') + theme +
  labs(y = 'Number of rhyGenes', x = 'Number of randomly selected samples')
  
p2 <- ggplot(data = a, aes(x = random.sample.size, y = rhyQTL.num)) + 
  geom_point() +
  geom_line(color = 'blue') + theme +
  labs(y = 'Number of rhyQTLs', x = 'Number of randomly selected samples')

pp1 <- ggplot(data = b, aes(x = random.sample.size, y = rhyGene.num)) + 
  geom_point() +
  geom_line(color = 'blue') + theme +
  labs(y = 'Number of rhyGenes', x = 'Number of randomly selected samples')

pp2 <- ggplot(data = b, aes(x = random.sample.size, y = rhyQTL.num)) + 
  geom_point() +
  geom_line(color = 'blue') + theme +
  labs(y = 'Number of rhyQTLs', x = 'Number of randomly selected samples')
ppp1 <- ggplot(data = c, aes(x = random.sample.size, y = rhyGene.num)) + 
  geom_point() +
  geom_line(color = 'blue') + theme +
  labs(y = 'Number of rhyGenes', x = 'Number of randomly selected samples')

ppp2 <- ggplot(data = c, aes(x = random.sample.size, y = rhyQTL.num)) + 
  geom_point() +
  geom_line(color = 'blue') + theme +
  labs(y = 'Number of rhyQTLs', x = 'Number of randomly selected samples')
pdf('../13_Figures/11_downsample.pdf', width = 5, height = 10, useDingbats = F)
print(p2, vp=viewport(.42, .18, x = .3, y = .8))
print(p1, vp=viewport(.4, .18, x = .7, y = .8))
print(pp2, vp=viewport(.42, .18, x = .3, y = .5))
print(pp1, vp=viewport(.4, .18, x = .7, y = .5))
print(ppp2, vp=viewport(.42, .18, x = .3, y = .3))
print(ppp1, vp=viewport(.4, .18, x = .7, y = .3))
dev.off()
