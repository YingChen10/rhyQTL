setwd('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/00_data/CPM_covariate_remove')
expression <- read.table('Adipose-Visceral_Omentum.txt', header = T)

random.size <- c(100, 150, 200, 250, 300, 350, 400, 450)

for(size in random.size){
  
  message(size)
  random.sample <- sample(names(expression), size)
  tmp <- select(expression, random.sample)  
  
  write.table(tmp, paste0('Adipose-Visceral_Omentum_subsample.', size, '.txt'), quote = F, sep = '\t')
  
}



