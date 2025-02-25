
donor.id <- fread('/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/head_sub.txt') %>% as.data.frame()


donor.id[,10:847]
shuffled.id <- as.character(c(donor.id[,1:9], sample(donor.id[,10:847])))

write.table(shuffled.id, '/workspace/rsrch1/ychen/Projects/Project03_human_circadian/rQTL/head_sub_shuffle1.txt', quote = F, sep ='\t', row.names = F) 
