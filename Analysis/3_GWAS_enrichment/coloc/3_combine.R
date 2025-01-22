
suppressMessages(library(tidyverse))
suppressMessages(library(data.table))
suppressMessages(library(parallel))
suppressMessages(library(gridExtra))
suppressMessages(library(grid))

args <- commandArgs(trailingOnly = TRUE)
tissue <- args[1]
tissue <- 'Liver'

dir <- "/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/"
cat <- "cis_QTL/"

################################################## 1. combine all pp4 0.4
indir <- paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/coloc/01_rhyQTL_coloc_PP4_filter/', tissue, '/')
file <- list.files(indir)

all <- c()
for(x in file){

  tmp <- read.table(paste0(indir, x), header = T)
  tmp$gwas.min.p <- as.numeric(tmp$gwas.min.p)
  tmp$PP4 <- as.numeric(tmp$PP4)
  all <- rbind(all, tmp)

}


all <- filter(all, gwas.min.p < 1e-8 & PP4 > 0.7)

indir <- paste0(dir, cat, 'Result/03_enrich_GWAS/coloc/')
# gene position annotation used to annotation gene in the plot
annot <- read.table(paste0(indir, 'gene.annot.txt'), header = T)

all <- merge(all, annot, by.x = 'rGene', by.y = 'gene_id')

filter <- all

names(filter)[names(filter) == "rGene"] <- 'gene'
names(filter)[names(filter)=='file'] <- 'split_pos_file'
names(filter)[names(filter)=='GWAS'] <- 'gwas'
filter$gwas <- gsub('EBI', 'ebi', filter$gwas)

gwas.list <- fread(paste0(dir, cat, './03_enrich_GWAS/GWAS_enrichment/GWAS_r2024-01-19.filter.txt'))
trait <- unique(select(gwas.list, Study_accession, Trait_mod))
names(trait) <- c('gwas', 'trait')
gwas.list2 <- fread(paste0(dir, cat, './Result/03_enrich_GWAS/coloc/GWAS.id.list.all'))
trait2 <- unique(select(gwas.list2, id, trait))
names(trait2) <- c('gwas', 'trait')
trait.list <- rbind(trait, trait2)

# to.plot <- merge(filter, trait.list, by = 'gwas') 
to.plot <- merge(filter, trait.list, by = 'gwas') 
to.plot$genome <- '38'


outdir <- paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/coloc/04_rhyQTL_PP4_filter_plot/')
if (!file.exists(outdir)) {dir.create(outdir, recursive = T)}

write.table(to.plot, paste0(outdir, tissue, '_04_PP4_filter_plot.txt'), quote = F, sep = '\t', row.names = F)

message('Combining PP4 results is done!')



######################################################### 2. plot 
Temp.dir <- paste0(dir, cat, 'Result/03_enrich_GWAS/coloc/Temp/')
dir.create(Temp.dir, showWarnings = FALSE, recursive = TRUE)  
Sys.setenv(TMPDIR = Temp.dir) 


source('~/Projects/Project03_human_circadian/Scripts/rQTL/cis_rQTL/Script/03_enrich_GWAS/2024-01/20241017_coloc/Function/coloc_plot.R')

tissue_indir <- tissue
indir <- paste0(dir, cat, 'Result/03_enrich_GWAS/coloc/')


file_name <- paste0(outdir, tissue, '_04_PP4_filter_plot.txt')
to.plot <- fread(file_name)


# Import TSS pos this is used for get eQTLs within the TSS flanking 1Mb region
tss <- read.table(paste0(indir, "Gene.TSS.GRCh38.txt"), head = T, sep="\t")
tss$Chr <- as.numeric(tss$Chr)

# match the name of rQTL and eQTL files #
name_match <- fread(paste0(indir, "file_name.match.txt"), head=T)
rQTL_tissue <- name_match$rQTL_file[name_match$rQTL_file == tissue_indir]
eQTL_tissue <- name_match$eQTL_file[name_match$rQTL_file == tissue_indir]

# read maf and rsid information
maf <-fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Data/03_MAF059.VCF.MAF.rsid.txt')
names(maf)[names(maf)=='rsid'] <- 'snp'

# gene position annotation used to annotation gene in the plot
annot <- read.table(paste0(indir, 'gene.annot.txt'), header = T)
outdir <- paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/coloc/04_rhyQTL_PP4_filter_plot/', tissue, '/')
plotoutdir <- paste0('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/cis_QTL/Result/03_enrich_GWAS/coloc/04_rhyQTL_PP4_filter_plot/plot/', tissue, '/')


if (!file.exists(outdir)) {dir.create(outdir)}
if (!file.exists(plotoutdir)) {dir.create(plotoutdir)}

lapply(1:nrow(to.plot), function(i){
  
  message(paste0('starting ', i))
  #print(to.plot[i,])
  
  # get information of gene tissue and rQTL file 
  genes <- to.plot$gene[i]
  split_pos_file <- to.plot$split_pos_file[i]
  #gwas_id <- to.plot$gwas[i]
  #message(paste0('processing ', i))
  gwas_id <- to.plot$gwas[i]
  trait <- to.plot$trait[i]
  name <- to.plot$name[i]
  
  file.name <- paste0(plotoutdir, paste(rQTL_tissue, gwas_id, gsub(' ', '_', trait), genes, name, 'png', sep = '.'))
  if (file.exists(file.name)) {
    message("exists")
    return(NULL)
  }
  
  # Import and format rQTL data
  rQTL_file <-  paste0(dir, cat, "00_rQTL_mapping/04_hanova/", rQTL_tissue, '/', split_pos_file, ".rds") 
  rqtl_list <- readRDS(rQTL_file)
  
  gene.tmp <- rqtl_list[[genes]][,-7]           
 
  gene.tmp$HANOVA.norm <- as.numeric(gene.tmp$HANOVA.norm)
  rqtl <- merge(gene.tmp, maf, by = 'ID')  
  duplicated_rows <- duplicated(rqtl$snp)
  rqtl <- rqtl[!duplicated_rows, ]
  rqtl$chr=as.numeric(gsub("chr","",rqtl$chr))
  rqtl=arrange(rqtl,pos)
  names(rqtl)[names(rqtl) == 'HANOVA.norm'] <- 'P.value'
  
  # get tss position
  tss_tmp <- filter(tss, Gene == genes)
  names(tss_tmp) <- c('gene', 'chr', 'pos')
  
  # Import and format eQTL data
  eQTL_file <- paste0(dir, cat, "03_enrich_GWAS/gwas/gtex.v8.eQTL/", eQTL_tissue, ".allpairs.txt.gz")
  eqtl <- Import_Format_eQTL(eQTL_file, genes, maf, tss_tmp)
  
  
  genome <- to.plot$genome[i]
  # convert tss position of the gene to hg19 version as GWAS information is based on hg19 version
  if(genome == '37'){
    tss_hg19 <- liftover_hg38_to_hg19(tss_tmp)
    tss_tmp$pos <- tss_hg19$pos
    #region <- paste0(tss_hg19$chr, ":", tss_hg19$pos - 1000000, "-", tss_hg19$pos + 1000000)
  }
  
  
  gwas_file <- paste0(indir ,"00_Data_GWAS_input/", gwas_id, ".sort.txt.gz")
  gwas <- fread(gwas_file, head=T)
  gwas$chr <- as.numeric(gwas$chr)
  gwas$pvalues <- as.numeric(gwas$pvalues)
  # remove duplicated snp
  duplicated_rows <- duplicated(gwas$snp)
  gwas <- gwas[!duplicated_rows, ]
  gwas <- na.omit(gwas)
  gwas <- Import.Format.GWAS(gwas, tss_tmp)
  
  if(genome == '37'){gwas <- liftover_hg19_to_hg38(gwas)}
  
  # get position and P value (take log) to plot
  rQTL_plot <- transmute(rqtl, pos, log10P=-log10(P.value))
  eQTL_plot <- transmute(eqtl, pos, log10P=-log10(P.value))
  gwas_plot <- transmute(gwas, pos, log10P=-log10(pvalues))
  
  
  rQTL_plot$pos <- as.numeric(rQTL_plot$pos)
  rQTL_plot$log10P <- as.numeric(rQTL_plot$log10P)
  
  eQTL_plot$pos <- as.numeric(eQTL_plot$pos)
  eQTL_plot$log10P <- as.numeric(eQTL_plot$log10P)
  
  gwas_plot$pos <- as.numeric(gwas_plot$pos)
  gwas_plot$log10P <- as.numeric(gwas_plot$log10P)
  
  ymax <- max(max(eQTL_plot$log10P), max(rQTL_plot$log10P))
  
  # gene annotation infor
  annot_tmp <- Gene.Annot(annot, tss_tmp, ymax)
  
  
  coloc = c(gwas_id, trait, genes, name, split_pos_file, to.plot$PP4[i], ymax)
  list <- list(gwas = gwas_plot, eqtl = eQTL_plot, rqtl = rQTL_plot, annot = annot_tmp, coloc = coloc)
  save(list, file = paste0(outdir, paste(rQTL_tissue, gwas_id, gsub(' ', '_', trait), genes, name, "RData", sep = '.')))
  
  png(paste0(plotoutdir, paste(rQTL_tissue, gwas_id, gsub(' ', '_', trait), genes, name, 'png',sep = '.')), width = 1800, height = 1800, res = 300)
  PLOT(gwas_plot, tss_tmp, eQTL_plot, ymax, rQTL_plot, annot_tmp)  
  grid.text(paste0(gwas_id, ' ', trait), x = unit(0.5, "npc"), y = unit(0.98, "npc"), gp = gpar(fontsize = 10))
  grid.text(paste0(rQTL_tissue, ' ', genes, ' ', name), x = unit(0.5, "npc"), y = unit(0.95, "npc"), gp = gpar(fontsize = 10))
  grid.text(signif(to.plot$PP4[i], digits = 3), x = unit(0.5, "npc"), y = unit(0.92, "npc"), gp = gpar(fontsize = 10))
  dev.off()
  
  message(paste0('ending ', i))
  
})


