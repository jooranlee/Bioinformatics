#summarizing, analyzing and visualizing MAF files
library(maftools)
library(R.utils)
library(readxl)
library(ggplot2)
library(dplyr)
#set working directory to where the files are located
setwd('/Users/jooranlee/Desktop/snu_BI/maftools')

#add data info
maf.somatic = read.maf('OVCA113_somatic.maf.gz')
loh <- read_excel('OVCA_WES.xlsx') 
loh$batch <- 'batch_1'
loh$batch[ loh$ID %in% sprintf("OVCA%03d", seq(89,113)) ] <- "batch_2"
colnames(loh) <- c('ID', 'purity', 'LOH', 'batch')

#boxplot
ggplot(loh) + geom_boxplot(aes(x=batch, y=LOH), outlier.shape = NA, notch = T) + 
  geom_jitter(aes(x=batch, y = LOH), width = 0.2)

#merge LOH with MAF summary 
tmp <- mafSummary(maf.somatic)

#summarizing MAF
getSampleSummary(maf.somatic)
getGeneSummary(maf.somatic)

plotmafSummary(maf.somatic)

#oncoplots
aml_genes <- c('TP53', 'EGFR', 'BRCA1', 'BRCA2','PTPN23', 'SLC25A5')
oncoplot(maf.somatic, genes = aml_genes) 
#group by pathways - oncoplots
oncoplot(maf.somatic, gene_mar = 8, fontSize = 1)

#transitions and trasversions 
titv.somatic <- titv(maf.somatic, plot = F, useSyn = T)
plotTiTv(titv.somatic)

#two samples are showing different ti/tv ratio
filter(titv.somatic$TiTv.fractions, Ti > 50)

#samples with exceptional C>T conversions: 
filter(titv.somatic$fraction.contribution, `C>T` > 30)

#labelling plots for selected genes
lollipopPlot(maf.somatic, gene = 'BRCA2', AACol = 'Protein_Change', showMutationRate = T)
lollipopPlot(maf.somatic, gene = 'BRCA1', AACol = 'Protein_Change', showMutationRate = T, labelPos = 881)

#rainfall plots
rainfallPlot(maf.somatic, detectChangePoints = TRUE, pointSize = 0.4)

#compare mutation load aginst TCGA cohorts
laml.mutload = tcgaCompare(maf.somatic, cohortName = 'marklogen', logscale = TRUE, 
                           capture_size = 50) 

#somatic interactions
#exclusive/co-occurence event analysis on top 10 mutated genes
somaticInteractions(maf.somatic, top = 25, pvalue = c(0.05, 0.1))

#detecting cancer driver genes based on positional clustering
laml.sig = oncodrive(maf.somatic, minMut = 5, 
                     pvalMethod = 'zscore')
head(laml.sig)
plotOncodrive(res = laml.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

#plotVaf
plotVaf(maf.somatic, vafCol = 'tumor_f')
