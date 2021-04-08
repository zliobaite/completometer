# 2020 10 31 I.Zliobaite
# find unique species of the Neogene

data_pb <- read.csv('pbdb_20210216_Miocene.csv', header = TRUE, sep = ",")
print('PBDB')
print(dim(data_pb))
ind <- which(data_pb[,'accepted_rank'] == 'species')
data_pb <- data_pb[ind,]
un_pb <- unique(data_pb[,'accepted_name'])
print(length(un_pb))
gen_pb <- c()
sp_pb <- c()
synonyms_pb <- c()
for (sk in 1:length(un_pb)){
  gen_now <- strsplit(as.vector(un_pb)[sk]," ")[[1]][1]
  gen_pb <- c(gen_pb,gen_now)
  nn <- length(strsplit(as.vector(un_pb)[sk]," ")[[1]])
  sp_now <- strsplit(as.vector(un_pb)[sk]," ")[[1]][nn]
  sp_pb <- c(sp_pb,sp_now)
  if (nn==3){
    syn_now <- strsplit(as.vector(un_pb)[sk]," ")[[1]][2]
    syn_now <- gsub("\\(|\\)", "", syn_now)
    synonyms_pb <- c(synonyms_pb,paste(syn_now,sp_now))
  }
}

un_sp_pb <- paste(gen_pb,sp_pb)
print(length(un_sp_pb))


data_all <- read.csv('NOW_20210216_public.csv', header = TRUE, sep = "\t")
print('NOW')
#ind <- which(data_all[,'LOC_STATUS'] =='public')
#data_all <- data_all[ind,]
#print(dim(data_all))

ind <- which(data_all[,'MIN_AGE']<23.3)
data_all <- data_all[ind,]
ind <- which(data_all[,'MAX_AGE']>5.333)
data_all <- data_all[ind,]

print(dim(data_all))

sum_sp <- 0

#ind_indet <- which(data_all[,'SPECIES']=='indet.')
ind <- which(data_all[,'SPECIES']!='indet.')
data_all <- data_all[ind,]
ind <- which(data_all[,'SPECIES']!='Indet.')
data_all <- data_all[ind,]
ind <- which(data_all[,'SPECIES']!='indet')
data_all <- data_all[ind,]

ind <- which(data_all[,'SPECIES']=='sp.')
sum_sp <- sum_sp + length(ind)
ind <- which(data_all[,'SPECIES']!='sp.')
data_all <- data_all[ind,]
ind <- which(data_all[,'SPECIES']=='sp')
sum_sp <- sum_sp + length(ind)
ind <- which(data_all[,'SPECIES']!='sp')
data_all <- data_all[ind,]
ind <- which(data_all[,'SPECIES']=='sp.')
sum_sp <- sum_sp + length(ind)
ind <- which(data_all[,'SPECIES']!='sp.')
data_all <- data_all[ind,]
ind <- which(data_all[,'SPECIES']=='sp. 1')
sum_sp <- sum_sp + length(ind)
ind <- which(data_all[,'SPECIES']!='sp. 1')
data_all <- data_all[ind,]
ind <- which(data_all[,'SPECIES']=='sp. a')
sum_sp <- sum_sp + length(ind)
ind <- which(data_all[,'SPECIES']!='sp. a')
data_all <- data_all[ind,]
ind <- which(data_all[,'SPECIES']=='sp. b')
sum_sp <- sum_sp + length(ind)
ind <- which(data_all[,'SPECIES']!='sp. b')
data_all <- data_all[ind,]
ind <- which(data_all[,'SPECIES']=='sp. c')
sum_sp <- sum_sp + length(ind)
ind <- which(data_all[,'SPECIES']!='sp. c')
data_all <- data_all[ind,]

ind <- which(data_all[,'GENUS']=='gen.')
sum_sp <- sum_sp + length(ind)
ind <- which(data_all[,'GENUS']!='gen.')
data_all <- data_all[ind,]
ind <- which(data_all[,'GENUS']=='Gen.')
sum_sp <- sum_sp + length(ind)
ind <- which(data_all[,'GENUS']!='Gen.')
data_all <- data_all[ind,]

ind <- which(data_all[,'GENUS']!='indet.')
data_all <- data_all[ind,]
ind <- which(data_all[,'GENUS']!='Indet.')
data_all <- data_all[ind,]

ind <- which(data_all[,'GENUS']=='incertae sedis')
sum_sp <- sum_sp + length(ind)
ind <- which(data_all[,'GENUS']!='incertae sedis')
data_all <- data_all[ind,]

un_sp_id <- as.vector(unique(data_all[,'SIDNUM']))
sp_names_nowdb <- paste(data_all[,'GENUS'],data_all[,'SPECIES'])
un_sp_nowdb <- unique(sp_names_nowdb)
print(length(un_sp_nowdb))

synonyms_nowdb <- c()
syn <- data_all[,'SYNONYMS']
for (sk in 1:length(syn)){
  syn_now <- syn[sk]
  str_now <- strsplit(as.vector(syn_now),":")[[1]]
  ll <- length(str_now)
  for (sk2 in 1:ll){
    ss_now <- trimws(str_now[sk2])
    ss_parts <- strsplit(as.vector(ss_now)," ")[[1]]
    if (length(ss_parts)>1){
      ss_sp_now <- paste(ss_parts[1],ss_parts[2])
      synonyms_nowdb <- c(synonyms_nowdb,ss_sp_now)
    }
  }
}


print('species at both')
sp_both <- intersect(un_sp_nowdb,un_sp_pb)
print(length(sp_both))
print('species in PBDB but not in NOW')
sp_only_pb <- setdiff(un_sp_pb,un_sp_nowdb)
print(length(sp_only_pb))
print('species in PBDB but not in NOW excluding NOW synonyms')
sp_only_pb <- setdiff(sp_only_pb,synonyms_nowdb)
print(length(sp_only_pb))
write.table(sp_only_pb, file = "out_species_only_pb.csv",sep = '\t')
print('species in PBDB but not in NOW excluding NOW and PBDB synonyms')
nsyn <- length(intersect(synonyms_pb,un_sp_nowdb))
print(length(sp_only_pb)-nsyn)


data_un <- c()
for (sk in 1:length(un_sp_nowdb)){
  ind <- which(sp_names_nowdb==un_sp_nowdb[sk])
  data_un <- rbind(data_un,data_all[ind[1],33:87])
}
write.table(data_un, file = "out_species_species.csv",sep = '\t')


data_tod <- read.csv('mdd_20200924.csv', header = TRUE, sep = ",")
ind <- which(data_tod[,'extinct']==0)
data_tod <- data_tod[ind,]
ind <- which(data_tod[,'domestic']==0)
data_tod <- data_tod[ind,]
library(stringr)
orders_living_all <- str_to_title(as.vector(data_tod[,'order']))
un_orders_living <- unique(orders_living_all)

un_orders_nowdb <- as.vector(unique(data_un[,'ORDER']))

un_orders_joint <- union(un_orders_living,un_orders_nowdb)

completometer <- c()
for (sk in 1:length(un_orders_joint)){
  order_now <- un_orders_joint[sk]
  count_nowdb <- length(which(data_un[,'ORDER']==order_now))
  count_living <- length(which(orders_living_all==order_now))
  completometer <- rbind(completometer,c(order_now,count_living,count_nowdb))
}

colnames(completometer) <- c('Order','no. living','no NOW')

write.table(completometer, file = "out_completometer.csv",sep = '\t')
