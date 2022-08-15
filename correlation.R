#see the correlation between PDC(patient-derived cells) and parental tissue using their WTS data, respectively
setwd("/Users/jooranlee/Dropbox/Mac/Desktop/correlation") 
df <- read.csv('dds_log2_tranform_220715.csv')
rownames(df) <- df$X
df$X <- NULL
#sample 4 and 6 were removed after PCA 
sample_no <- c('sample2_p', 'sample3_p', 'sample5_p', 'sample7_p', 'sample1_p', 'sample2_t', 'sample1_t', 'sample3_t', 'sample5_t', 'sample7_t')
colnames(df) <- sample_no

library(ggplot2)
library(reshape2)
library(dplyr)
df<-df[c(paste0(paste0("sample", c(1,2,3,5,7)), "_p"), paste0(paste0("sample", c(1,2,3,5,7)), "_t"))]
sample.cor <- cor(df, method = c('spearman'))
gene.cor<-apply(df, MARGIN=1, FUN=function(x){return(cor(as.numeric(x[1:5]), as.numeric(x[6:10]), method = "spearman") )}) 
#no filtering (30867 genes)
ggplot(melt(sample.cor) %>% dplyr::filter(melt(sample.cor)$Var1 %in% grep(sample_no, pattern="_t$", value=T), 
                     melt(sample.cor)$Var2 %in% grep(sample_no, pattern="_p$", value=T))) + geom_tile(aes(x=Var1, y=Var2, fill=value))
#check the histogram
hist(melt(gene.cor)$value)

#0.9 threshold (1437 genes)
sample.cor<-cor(df[gene.cor>0.9, ], method="spearman")
ggplot(melt(sample.cor) %>% dplyr::filter(melt(sample.cor)$Var1 %in% grep(sample_no, pattern="_t$", value=T), 
                                          melt(sample.cor)$Var2 %in% grep(sample_no, pattern="_p$", value=T))) + geom_tile(aes(x=Var1, y=Var2, fill=value))

#0.7 threshold (7350 genes) 
sample.cor<-cor(df[gene.cor>0.7, ], method="spearman")
ggplot(melt(sample.cor) %>% dplyr::filter(melt(sample.cor)$Var1 %in% grep(sample_no, pattern="_t$", value=T), 
              melt(sample.cor)$Var2 %in% grep(sample_no, pattern="_p$", value=T))) + geom_tile(aes(x=Var1, y=Var2, fill=value))


#0.5 threshold (14346 genes)
sample.cor<-cor(df[gene.cor>0.5, ], method="spearman")
ggplot(melt(sample.cor) %>% dplyr::filter(melt(sample.cor)$Var1 %in% grep(sample_no, pattern="_t$", value=T), 
              Var2 %in% grep(sample_no, pattern="_p$", value=T))) + geom_tile(aes(x=Var1, y=Var2, fill=value))


#0.3 threshold (18093 genes)
sample.cor<-cor(df[gene.cor>0.3, ], method="spearman")
ggplot(melt(sample.cor) %>% dplyr::filter(melt(sample.cor)$Var1 %in% grep(sample_no, pattern="_t$", value=T), 
              Var2 %in% grep(sample_no, pattern="_p$", value=T))) + geom_tile(aes(x=Var1, y=Var2, fill=value))

#0.1 threshold (22295 genes)
sample.cor<-cor(df[gene.cor>0.1, ], method="spearman")
ggplot(melt(sample.cor) %>% dplyr::filter(melt(sample.cor)$Var1 %in% grep(sample_no, pattern="_t$", value=T), 
              Var2 %in% grep(sample_no, pattern="_p$", value=T))) + geom_tile(aes(x=Var1, y=Var2, fill=value))




