setwd('F:/Projects/Project03_human_circadian/rQTL/cis_QTL/')
expression <- fread('../GTEx_nor_expression/Adipose-Visceral_Omentum.txt')
# import expression file
expression <- as.data.frame(expression)
# expression <- expression[,2:(ncol(expression) - 34)]
# # expression$EnsemblID <- sapply(strsplit(expression$EnsemblID, '_'), "[", 1)
# # names(expression)[2:ncol(expression)] <- sapply(strsplit(names(expression)[2:ncol(expression)], '\\.'), "[", 2)
# # 
# # head <- read.table('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/head_sub.txt')
# 
# expression <- select(expression, EnsemblID, intersect(as.character(head[1,]), names(expression)))
# 

random.size <- c(100, 150, 200, 250, 300, 350, 400, 450)

for(size in random.size){
  
  message(size)
  random.sample <- sample(names(expression[,2:(ncol(expression) - 34)]), size)
  tmp <- select(expression, EnsemblID, random.sample)  
  tmp <- cbind(tmp, expression[,c((ncol(expression) - 33) : ncol(expression))])
  write.table(tmp, paste0('../GTEx_nor_expression/Adipose-Visceral_Omentum_subsample.', size, '.txt'), quote = F, sep = '\t', row.names = F)
  
}


