# We obtained the single-cell transcriptome data about Gata1+ states of hematopoietic progenitor cells (HPCs) from Wang et al. (https://cospar.readthedocs.io/).
# We first excluded the cells without clonal information and changed the dataset into the CSV format according to different types of annotation.

# Take the data of day 2 as an example of LUG searching algorithm

rm(list=ls())

clone_inf <- read.csv(file='~/.../clone_Hematopoiesis.csv',row.names=1,header=T)
mat <- read.csv(file='~/.../Hematopoiesis_data.csv',row.names=1,header=T)
time_inf <- read.csv(file='~/.../Hematopoiesis_obs.csv',row.names=1,header=T)
biasgenes_inf <- read.csv(file='~/.../bias_gene.csv',header=T)
biasgenes_vector_1 <- c(unlist(biasgenes_inf[,1]),unlist(biasgenes_inf[,2]))
biasgenes_vector_1 <- setdiff(biasgenes_vector_1,"")
names(biasgenes_vector_1) <- c(rep(colnames(biasgenes_inf)[1],length(setdiff(unlist(biasgenes_inf[,1]),''))),rep(colnames(biasgenes_inf)[2],length(unlist(biasgenes_inf[,2]))))

mat <- t(as.matrix(mat))
mat <- as.data.frame(log1p(t(t(mat) / colSums(mat)) * 10000))

# Select the data of day 2
mat_d2 <- mat[,intersect(rownames(time_inf[time_inf$time_info=='2',]),colnames(mat))]
clone_inf_d2 <- clone_inf[colnames(mat_d2),]

biasgenes <- rowSums(as.matrix(mat_d2))
biasgenes <- biasgenes[biasgenes!=0]
biasgenes_vector <- names(biasgenes)

# Candidates of LUGs on day 2
Np_biasgene <- vector()
cellnum_biasgene <- vector()

for (n in 1:length(biasgenes_vector)) {
  
  zzdata <- biasgenes_vector[n]
  idata <- as.vector(unlist(mat_d2[zzdata,]))
  names(idata) <- colnames(mat_d2)
  idata <- sort(idata,decreasing=T)
  idata <- idata[idata!=0]
  idata_names <- names(idata)
  
  target_cell <- vector()
  
  if (length(idata)<ceiling(0.1*ncol(mat_d2))) {
    
    target_cell <- append(target_cell,idata_names)
    
  } else {
    
    idata_names <- idata_names[1:ceiling(0.1*ncol(mat_d2))]
    target_cell <- append(target_cell,idata_names)
    
  }
  
  cellnum_biasgene <- append(cellnum_biasgene,length(target_cell))
  
  clone_inf_target <- clone_inf_d2[target_cell,]
  clone_inf_targetlist <- vector('list',ncol(clone_inf_target))
  names(clone_inf_targetlist) <- colnames(clone_inf_target)
  
  for (i in 1:ncol(clone_inf_target)) {
    
    idata_i <- clone_inf_target[,i]
    names(idata_i) <- rownames(clone_inf_target)
    clone_inf_targetlist[i] <- list(names(idata_i[which(idata_i=='True')]))
    
  }
  
  cr <- lengths(clone_inf_targetlist)[lengths(clone_inf_targetlist)!=0]
  cr_event <- choose(cr,2)
  tmpdata <- sum(cr_event)/choose(length(target_cell),2)
  Np_biasgene <- append(Np_biasgene,1/tmpdata)
  
}

df <- data.frame(biasgenes_vector,Np_biasgene,cellnum_biasgene)
names(df) <- c('biasgenes','Np','cellnum')

# 1,000 times of ramdomly sampling 10% of cells for LUG searching algorithm 
Np_sh <- vector()

for (n in 1:1000) {
  
  target_cell <- sample(colnames(mat_d2),ceiling(0.1*ncol(mat_d2)),prob=rep(1,ncol(mat_d2)),replace=FALSE)
  clone_inf_target <- clone_inf_d2[target_cell,]
  clone_inf_targetlist <- vector('list',ncol(clone_inf_target))
  names(clone_inf_targetlist) <- colnames(clone_inf_target)
  
  for (i in 1:ncol(clone_inf_target)) {
    
    idata_i <- clone_inf_target[,i]
    names(idata_i) <- rownames(clone_inf_target)
    clone_inf_targetlist[i] <- list(names(idata_i[which(idata_i=='True')]))
    
  }
  
  cr <- lengths(clone_inf_targetlist)[lengths(clone_inf_targetlist)!=0]
  cr_event <- choose(cr,2)
  tmpdata <- sum(cr_event)/choose(length(target_cell),2)
  Np_sh <- append(Np_sh,1/tmpdata)
  
}

Np_sh <- sort(Np_sh)
Np_sh[1:10]
# [1] 598.9286 645.0000 698.7500 698.7500 762.2727 838.5000 838.5000 838.5000
# [9] 838.5000 838.5000

# P-value < 0.01 and LUG is expressed in at least 10 percent of the cells 
df_1 <- df[is.finite(df$Np),]
df_1 <- df_1[df_1$Np<838.5,]
df_1 <- df_1[df_1$cellnum>64,]

# LUGs of day 2
df_1$biasgenes