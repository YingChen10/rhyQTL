
library(parallel)
dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/"
cat <- "cis_QTL/"

# file <- list.files(paste0(dir, 'GTEx_nor_expression/'), pattern = 'txt')
# tissue <- gsub('.txt', '', file)
# tissue.list <- tissue[tissue != 'Kidney-Cortex']
# tissue.list <- tissue.list[!grepl('female', tissue.list) & !grepl('male', tissue.list)]
file.list <- list.files(paste0(dir, cat, 'split_pos'))


args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
indir <- paste0(dir, cat, '00_rQTL_mapping/00_Genotype/')

all <- mclapply(file.list, function(file){
  
  
  tmp <- readRDS(paste0(indir, tissue, '/', file, '.rds'))
  pval <- lapply(1:length(tmp), function(i){
    
    gene <- names(tmp)[i]
    gene.tmp <- tmp[[i]]
    
    
    return(data.frame(
      SNP = names(gene.tmp),
      gene = gene,
      pvalue = 1)) # Here we assign 1 to each locus in order to using eigenMT to get Meff
    
  })
  
  pval <- as.data.frame(do.call(rbind, pval)) 
  
}, mc.cores = 5)

all <- as.data.frame(do.call(rbind, all))

write.table(all, paste0(dir, cat, '00_rQTL_mapping/00_Meff_input/', tissue, '.txt'), quote = F, sep = '\t', row.names = F)

# ### generate input genotype and phenotype annotation file for Meff
# 
# ## 1. generate Meff_input files
# cd /workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL
# 
# grep -v "^##" GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hwe_MAF059.vcf | cut -f 1,10- 1 |sed "s/#//g" | sed -e 's/0|0/0/g' -e 's/0|1/1/g' -e 's/1|0/1/g' -e 's/1|1/2/g' -e 's/\./NA/g' > genotypes.txt
# 
# grep -v "^##" GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.SHAPEIT2_phased.MAF01.hwe_MAF059.vcf > 1
# sed '1d' 1 | cut -f 1,2,3 | awk 'BEGIN {print "snp\tchr_snp\tpos"} {print $3"\t"$1"\t"$2}' > gen.positions.txt
# 
# awk 'BEGIN {print "gene_id\tchrom_probe\ts1\ts2"}{if ($4 == "+") print $5"\t"$1"\t"$2"\t"$2;else if ($4 == "-") print $5"\t"$1"\t"$3"\t"$3}' all_gene.annot > phe.positions.txt
# 
# ## 2. split according to chromosomes
# sh genotypes.split.sh genotypes.txt genotypes
# sh genotypes.split.sh gen.positions.txt gen.positions 
# 
# # ==================== cat genotype.split.sh =======================
# #!/bin/bash
# header=$(head -n 1 $1)
# mkdir $2
# 
# tail -n +2 gen.positions.txt | cut -f 1 | awk -F'_' '{print $1}' | sort | uniq > unique_elements.txt
# 
# for elem in $(cat unique_elements.txt); do
# 
# echo  $elem
# echo -e "$header" > "./$2/${elem}.txt"
# grep -P "^$elem\_" $1 >> "./$2/${elem}.txt"
# done
# # =========================== end =================================


##################
#
# plot
#
#################
all <- data.frame()
file.list <- list.files('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/00_rQTL_mapping/00_Meff_output', pattern = '.txt')
for(x in file.list){
  
  
  meff <- fread(paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/00_rQTL_mapping/00_Meff_output/', x))
  meff$tissue <- gsub(".txt", "", x)
  
  all <- rbind(all, meff)
  
}

median <- as.data.frame(all %>% group_by(tissue) %>% summarise(median.meff = median(Meff), mean.meff = mean(Meff)))
median$p.cut.off <- 0.01/median$median.meff
median$p.cut.off.mean <- 0.01/median$mean.meff


theme <- theme_bw() +
  theme(axis.line = element_blank(),  
        axis.text = element_text(color = "black", size = 8),
        panel.grid = element_blank(),
        panel.border = element_rect(color = "black", fill = NA),
        axis.title = element_text(size = 8)) 

color <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Data/color_annotation_v3.txt')

all <- merge(all, color, by.x = 'tissue', by.y = 'tissue_rQTL')
sample.size <- read.table('Data/00_sample_size_with_phenotype.txt', header = T)
sample.size <- merge(sample.size, color, by.x = 'tissue', by.y = 'tissue_rQTL')
tissue.order <- sample.size[order(-sample.size$samplesize),]$FullName
all$FullName <- factor(all$FullName, levels = tissue.order)
write.table(all, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Results/00_Meff/meff.all.txt', quote = F, sep = '\t', row.names = F)


ggplot(data = all, aes(x = FullName, y = Meff)) + 
  geom_boxplot(outlier.shape = NA) + theme + ylim(0,40) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5)) +
  geom_hline(yintercept = median(all$Meff), linetype = 'dashed') 


all$Meff.2 <- all$Meff
all$Meff.2[all$Meff > 60] <- 60

p2 <- ggplot(data = filter(all, Meff < 60), aes(x = Meff)) + 
  geom_histogram(fill = 'skyblue', color = 'black', binwidth = 3) + theme + 
  geom_vline(xintercept = 15) 
p2

p1 <- ggplot(data = filter(all, Meff < 50), aes(x = Meff)) + 
  geom_density(fill = 'skyblue', color = 'black', adjust = 4) + 
  geom_vline(xintercept = 15) + theme

pdf("/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/13_Figures/00_Meff/meff_distribution_all.pdf",useDingbats = F, height = 2.4, width = 2.8)
print(p1)
print(p2)
dev.off()