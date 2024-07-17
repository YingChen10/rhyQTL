


file <- list.files('00_rQTL_mapping/04_Results/cis_rhyQTL_tissue/')
all <- c()

## define AM rhyGenes and PM rhyGenes
for(i in 1:length(file)){
  
  tissue <- gsub('.rQTL', '', file[i])
  
  tmp <- fread(paste0('00_rQTL_mapping/04_Results/cis_rhyQTL_tissue/', file[i]))
  tmp$phase <- apply(select(tmp, contains("phase"), max_amp_idx), 1, function(x) (x[x[3]]))
  mean.phase <- as.data.frame(tmp %>% group_by(gene) %>% summarise(mean.phase = median(phase)))
 
  
  annot <- fread('all_gene.annot')
  mean.phase <- merge(mean.phase, annot, by.x = 'gene', by.y = 'V6')
  mean.phase$tissue <- gsub('.rQTL', '', file[i])
  all <- rbind(all, mean.phase[,c(1,2,8,12)])
  
  write.table(filter(mean.phase, mean.phase > 0 & mean.phase <= 12)[,8], paste0('Result/00_rhyGene_phase/rGene_phase_12h_interval/', tissue, '.AM.rGene.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
  write.table(filter(mean.phase, mean.phase > 12 & mean.phase <= 24)[,8], paste0('Result/00_rhyGene_phase/rGene_phase_12h_interval/', tissue, '.PM.rGene.txt'), quote = F, sep = '\t', row.names = F, col.names = F)
  
}

write.table(all, 'Result/00_rhyGene_phase/all.tissue.rGene.phase.txt', quote = F, sep = '\t', row.names = F, col.names = F)  


###################################################### Plot all types of tissues ####################################
all <- fread('Result/00_rhyGene_phase/all.tissue.rGene.phase.txt')
names(all) <- c('gene', 'mean.phase', 'name', 'tissue')
color <- fread('Data/color_annotation_v3.txt')
all <- filter(all, tissue != 'Brain-Substantianigra')


for(x in unique(color$Class)[unique(color$Class) !=  "Cells"]){
  
  tmp <- filter(all, tissue %in% color$tissue_rQTL[color$Class == x])
  tmp$tissue <- factor(tmp$tissue, levels=unique(tmp$tissue))
  
  pdf(paste0('./Figure/00_rhyGene_phase/all_class.', gsub('/', '', x),'_phase_distribution.pdf'), width = 5, height = 2, useDingbats = F)
  print(
    ggplot(data = tmp, aes(x = mean.phase, color = tissue)) + geom_density() + coord_polar() + theme_minimal() +
      scale_colour_manual(values = color$ColorPlot[match(levels(tmp$tissue), color$tissue_rQTL)]) + 
      scale_x_continuous(limits = c(0, 24), breaks = c(0, 4, 8, 12, 16, 20, 24)) + 
      labs(x = 'Time of day (h)', y = 'Density')+
      theme(axis.line = element_blank(),
            axis.text = element_text(color = "black", size = 8),
            axis.title = element_text(size = 8),
            plot.title = element_text(size = 8)) +
      ggtitle(x)
    )
  dev.off()
}


###################################### Skin ############################################

m.1 <- fread('Result/00_rhyGene_phase/rGene_phase_12h_interval/Skin-SunExposed_Lowerleg.AM.rGene.txt', header = F)
m.2 <- fread('Result/00_rhyGene_phase/rGene_phase_12h_interval/Skin-NotSunExposed_Suprapubic.AM.rGene.txt', header = F)

m.common.am <- intersect(m.1$V1, m.2$V1)
sun.unique.am <- setdiff(m.1$V1, m.common.am)
no.sun.unique.am <- setdiff(m.2$V1, m.common.am)

write.table(m.common.am, 'Result/00_rhyGene_phase/skin/Skin-AM-common-rGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(sun.unique.am, 'Result/00_rhyGene_phase/skin/Skin-AM-SunExposed-unique-rGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(no.sun.unique.am, 'Result/00_rhyGene_phase/skin/Skin-AM-Nosun-unique-common-rGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)

m.1 <- fread('Result/00_rhyGene_phase/rGene_phase_12h_interval/Skin-SunExposed_Lowerleg.PM.rGene.txt', header = F)
m.2 <- fread('Result/00_rhyGene_phase/rGene_phase_12h_interval/Skin-NotSunExposed_Suprapubic.PM.rGene.txt', header = F)

m.common.pm <- intersect(m.1$V1, m.2$V1)
sun.unique.pm <- setdiff(m.1$V1, m.common.pm)
no.sun.unique.pm <- setdiff(m.2$V1, m.common.pm)

write.table(m.common.pm, 'Result/00_rhyGene_phase/skin/Skin-PM-common-rGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(sun.unique.pm, 'Result/00_rhyGene_phase/skin/Skin-PM-SunExposed-unique-rGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)
write.table(no.sun.unique.pm, 'Result/00_rhyGene_phase/skin/Skin-PM-Nosun-unique-common-rGene.txt', 
            sep = '\t', quote = F, row.names = F, col.names = F)



