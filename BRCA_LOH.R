library(DESeq2)
library(dplyr) 
library(readxl)
library(ggplot2)
getwd()
setwd("/Users/jooranlee/Desktop/BRCA/count")
x<- dir(path = ".", pattern = ".count$", recursive = T, full.names = T)
y<- read_excel('soamtic variants_113.xlsx')
sampleTable <- data.frame(sampleName = y$sample, directory = x)
sampleTable

#add type
type<- c(rep('type1', 88), rep('type2', 25)) 
sampleTable$type <- type
sampleTable$LOH <- y$LOH
sampleTable$BRCA_s_status <- y$BRCA_s_status
colnames(sampleTable) <- c('sample', 'directory', 'type', 'LOH', 'BRCA_s_status') 
#quantile
#25%: 15, 75%: 24
sampleTable<- mutate(sampleTable, loh_class = case_when(LOH>=24 ~ "hiLOH", LOH <= 15 ~"loLOH", TRUE ~ NA_character_))
sampleTable$type <- as.factor(sampleTable$type)
sampleTable$loh_class <- as.factor(sampleTable$loh_class)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(sampleTable = dplyr::filter(sampleTable, !is.na(loh_class)), directory = "/Users/jooranlee/Desktop/BRCA/count", 
                                       design = ~loh_class+type) 

assays(ddsHTSeq) #64

#plotPCA 
vsd2<-vst(ddsHTSeq)
plotPCA(vsd2, intgroup = 'type')
#input is a matrix of log transformed values
pca <- prcomp(t(assay(vsd2)))
tmp <- cbind(sampleTable, pca$x)
ggplot(tmp) + geom_point(aes(x= PC3, y = PC4, color = type))

##deseq
dds2 <- DESeq(ddsHTSeq)
summary(results(dds2, contrast = c("loh_class", "hiLOH", "loLOH")))
boxplot(log10(assays(dds2)[['cooks']]), intgroup = 'group')
    
dds.out <- results(dds2, contrast = c('loh_class', 'hiLOH', 'loLOH')) %>% as.data.frame()
dds.out %>% dplyr::filter(padj < 0.05) %>% dim()

###
library(tidyverse)
library(msigdbr)
library(clusterProfiler)
library(biomaRt)
organism = "org.Hs.eg.db"
library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(pathview)
library(AnnotationHub)
library(ensembldb)
library(AnnotationDbi)

dds.fgsea <- dds.out
ens <- rownames(dds.fgsea)
ensLookup <- gsub("\\.[0-9]*$", "", ens)
dds.fgsea$gene <- ensLookup
gene_id <- getBM(
  mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl")), 
  attributes = c("ensembl_gene_id", "gene_biotype", "external_gene_name", "entrezgene_id"), 
  filter = "ensembl_gene_id",
  values = ensLookup, 
  uniqueRows = TRUE)
colnames(gene_id) <- c('ensembl_gene_id', 'gene_biotype', 'external_gene_name', 'entrezgene_id')
dds2.fgsea <- left_join(dds.fgsea, gene_id[,c('ensembl_gene_id', 'gene_biotype', 'external_gene_name', 'entrezgene_id')], 
                        by = c('gene' = 'ensembl_gene_id'))

ens2symbol <- AnnotationDbi::select(org.Hs.eg.db, key = dds2.fgsea$gene, columns="SYMBOL", keytype="ENSEMBL")

ens2symbol <- as_tibble(ens2symbol)
dds2.fgsea <- inner_join(dds2.fgsea, ens2symbol, by = c('gene'='ENSEMBL'))
dds2.fgsea.sig <- dds2.fgsea %>% dplyr::filter(padj < 0.05) %>% 
  dplyr::select(SYMBOL, stat) %>% na.omit() %>% distinct()

#check how many overlapping genes are there
gene_list <- read_excel('gene_lists.xlsx') #230개의 HRD gene
overlapping_genes <- intersect(gene_list$`Gene Symbol`, dds2.fgsea.sig $SYMBOL) #3

library(fgsea)
library(tidyverse)
dds2.fgsea <- dds2.fgsea %>%
  dplyr::select(SYMBOL, stat) %>% 
  na.omit() %>%
  distinct()
ranks <- deframe(dds2.fgsea)
head(ranks, 20)

pathways.oncogene <- gmtPathways('/Users/jooranlee/Dropbox/Mac/Desktop/BRCA/count/c6.all.v7.5.1.symbols.gmt')
pathways.oncogene %>%
  head() %>% 
  lapply(head)
fgseaRes <- fgsea(pathways = pathways.oncogene, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>% 
  as_tibble() %>% 
  arrange(desc(NES))

ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) + 
  geom_col(aes(fill = padj<0.05)) + 
  coord_flip()+ 
  labs(x = 'pathway', y ='normalished ES', title = 'oncogenic pathway NES from GSEA')


###
tmp <- as.data.frame(assay(vsd2)) 
tmp$gid <- row.names(tmp)
tmp$gid<- gsub("\\.[0-9]*$", "", tmp$gid)
tmp<-left_join(tmp, ens2symbol, by = c("gid" = "ENSEMBL"))
library(ComplexHeatmap)
Heatmap(dplyr::filter(tmp, SYMBOL %in% gene_list$`Gene Symbol`)[,1:64])
