# ~ R script Gene Ontology analyses (Severo et al PNAS 2018 - scRNA-Seq of Anopheles gambiae hemocytes) ~ #

#### Load and prepare the data ####

## to load mapped read counts, prepare the data and subset cells 

counts <- read.table("~/Desktop/GitHub/20170126-count.header.star.tab", header=T,  row.names=1)

columns <- c("A1", "A2", "A4", "A5", "A6", "A7", "B8", "C1", "C3", "C7", "C8", "D1", "D3", "D4", "D5", "D9", "E4", "F1", "F3", "F5", "F6", "F7", "F9", "F10", "G1", "G3", "G4", "H4")
all_counts <- subset(counts, select=columns)

features <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")

all_counts_genes <- all_counts[-which(rownames(all_counts) %in% features), ]

all_counts_features <- all_counts[which(rownames(all_counts) %in% features), ]

## discard outlier cell (Figure S1) and subset single cells from pools

sc_counts <- subset(all_counts_genes, select=-c(A1, A2, D1))

## to subset endogenous and ERCC counts

geneTypes <- factor( c( AG="Ag", ER="ERCC" )[
  substr( rownames(sc_counts), 1, 2 ) ] )
countsAg <- sc_counts[ which( geneTypes=="Ag" ), ]
countsERCC <- sc_counts[ which( geneTypes=="ERCC" ), ]

## to divide by the size factors to get normalized counts (by sfERCC)

library("DESeq2")

sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsAg <- t( t(countsAg) / sfERCC)


#### Gene ontology (GO) analyses using topGO ####

library(topGO)
library(plyr)
library(genefilter)
source("https://bioconductor.org/biocLite.R")
biocLite("org.Ag.eg.db")

## GO terms in 90% of the cells (Table S4)

sc_nCountsAg1 <- nCountsAg[rowSums(nCountsAg) >=1,]
filter_sc_cells <- apply(sc_nCountsAg1, 1, function(x) length(x[x>=1])>=22)
all_sc_cells <- sc_nCountsAg1 [filter_sc_cells,]

## to run topGO - 90% of the cells

geneNames <- rownames(nCountsAg)
geneList <- factor(as.integer(geneNames %in% rownames(all_sc_cells)))
names(geneList) <- geneNames
str(geneList)

GOdata_all_sc_cells_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Ag.eg.db", ID ="ensembl")

resultTopGO.elim_all_sc_cells_BP <- runTest(GOdata_all_sc_cells_BP, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic_all_sc_cells_BP <- runTest(GOdata_all_sc_cells_BP, algorithm = "classic", statistic = "Fisher" )

topGOres.elim_all_sc_cells_BP <- GenTable(GOdata_all_sc_cells_BP, classic = resultTopGO.elim_all_sc_cells_BP , orderBy = "classic", topNodes=200)
topGOres.clas_all_sc_cells_BP <- GenTable(GOdata_all_sc_cells_BP, classic = resultTopGO.classic_all_sc_cells_BP , orderBy = "classic", topNodes=200)


#### Comparison PPO high and PPO low populations #####

samplenames2 <-c("Group 1", "Group 2",  "Group 1", "Group 1",	"Group 1",	"Group 2",	"Group 1", "Group 1",	"Group 2",	"Group 1", "Group 1", "Group 1", "Group 2",	"Group 1", "Group 1", "Group 1", "Group 1",	"Group 1", "Group 1", "Group 2", "Group 1", "Group 2", "Group 1")

## to filter out samples that are outliers on the PCA (Figure 2)

countsAg2 <- subset(countsAg, select=-c(A4, B8))

## to estimate size factors without the outlier samples

countsERCC2 <- subset(countsERCC, select=-c(A4, B8))
sfERCC2 <- estimateSizeFactorsForMatrix( countsERCC2 )

## to build the object and run the DESeq

conds = data.frame(samplenames2)
colnames(conds)="condition"
cds <- DESeqDataSetFromMatrix(countData = countsAg2, colData = conds, design = ~ condition)
sizeFactors(cds) <- sfERCC2
dds <- DESeq(cds, minReplicatesForReplace=Inf, fitType="local")
res=results( dds, contrast=c("condition","Group 2", "Group 1"), cooksCutoff=FALSE)
res<-res[order(res$padj),]

## to build a df with the DE genes and inspect expression of low padj ones

res.df <- as.data.frame(res)
Deseq.genes <- res.df[res.df$padj <= 0.1,]
Deseq.genes <-na.omit(Deseq.genes)

## make a data frame of the logFC results

Deseq.genes.pos <- Deseq.genes[Deseq.genes$log2FoldChange >= 0,]
Deseq.genes.neg <- Deseq.genes[Deseq.genes$log2FoldChange <= 0,]


## comparison of PPO populations based on GO terms (Table S6)

## to run topGO - Deseq pos

geneNames <- rownames(nCountsAg)
geneList <- factor(as.integer(geneNames %in% rownames(Deseq.genes.pos)))
names(geneList) <- geneNames
str(geneList)

GOdata_Deseq.genes.pos_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Ag.eg.db", ID ="ensembl")

resultTopGO.elim_Deseq.genes.pos_BP <- runTest(GOdata_Deseq.genes.pos_BP, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic_Deseq.genes.pos_BP <- runTest(GOdata_Deseq.genes.pos_BP, algorithm = "classic", statistic = "Fisher" )

topGOres.elim_Deseq.genes.pos_BP <- GenTable(GOdata_Deseq.genes.pos_BP, classic = resultTopGO.elim_Deseq.genes.pos_BP , orderBy = "classic", topNodes=200)
topGOres.clas_Deseq.genes.pos_BP <- GenTable(GOdata_Deseq.genes.pos_BP, classic = resultTopGO.classic_Deseq.genes.pos_BP , orderBy = "classic", topNodes=200)

## to run topGO - Deseq neg

geneNames <- rownames(nCountsAg)
geneList <- factor(as.integer(geneNames %in% rownames(Deseq.genes.neg)))
names(geneList) <- geneNames
str(geneList)

GOdata_Deseq.genes.neg_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Ag.eg.db", ID ="ensembl")

resultTopGO.elim_Deseq.genes.neg_BP <- runTest(GOdata_Deseq.genes.neg_BP, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic_Deseq.genes.neg_BP <- runTest(GOdata_Deseq.genes.neg_BP, algorithm = "classic", statistic = "Fisher" )

topGOres.elim_Deseq.genes.neg_BP <- GenTable(GOdata_Deseq.genes.neg_BP, classic = resultTopGO.elim_Deseq.genes.neg_BP , orderBy = "classic", topNodes=200)
topGOres.clas_Deseq.genes.neg_BP <- GenTable(GOdata_Deseq.genes.neg_BP, classic = resultTopGO.classic_Deseq.genes.neg_BP , orderBy = "classic", topNodes=200)


# sessionInfo()
# R version 3.3.2 (2016-10-31)
# Platform: x86_64-apple-darwin13.4.0 (64-bit)
# Running under: OS X El Capitan 10.11.6

# locale:
# [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

# attached base packages:
#  [1] stats4    parallel  stats     graphics  grDevices utils     datasets  methods   base     

# other attached packages:
# [1] topGO_2.26.0         SparseM_1.77         GO.db_3.4.0          AnnotationDbi_1.36.2
# [5] IRanges_2.8.2        S4Vectors_0.12.2     Biobase_2.34.0       graph_1.52.0        
# [9] BiocGenerics_0.20.0  biomaRt_2.30.0      

# loaded via a namespace (and not attached):
# [1] Rcpp_0.12.16               locfit_1.5-9.1             lattice_0.20-35           
# [4] digest_0.6.15              GenomeInfoDb_1.10.3        plyr_1.8.4                
# [7] backports_1.1.2            acepack_1.4.1              RSQLite_2.1.0             
# [10] ggplot2_2.2.1              pillar_1.2.1               zlibbioc_1.20.0           
# [13] rlang_0.2.0                lazyeval_0.2.1             rstudioapi_0.7            
# [16] data.table_1.10.4-3        annotate_1.52.1            blob_1.1.1                
# [19] rpart_4.1-13               Matrix_1.2-12              checkmate_1.8.5           
# [22] splines_3.3.2              BiocParallel_1.8.2         geneplotter_1.52.0        
# [25] stringr_1.3.0              foreign_0.8-69             htmlwidgets_1.0           
# [28] RCurl_1.95-4.10            bit_1.1-12                 munsell_0.4.3             
# [31] pkgconfig_2.0.1            base64enc_0.1-3            htmltools_0.3.6           
# [34] nnet_7.3-12                SummarizedExperiment_1.4.0 tibble_1.4.2              
# [37] gridExtra_2.3              htmlTable_1.11.2           matrixStats_0.53.1        
# [40] Hmisc_4.1-1                XML_3.98-1.10              bitops_1.0-6              
# [43] grid_3.3.2                 xtable_1.8-2               gtable_0.2.0              
# [46] DBI_0.8                    magrittr_1.5               scales_0.5.0              
# [49] stringi_1.1.7              XVector_0.14.1             genefilter_1.56.0         
# [52] latticeExtra_0.6-28        Formula_1.2-2              RColorBrewer_1.1-2        
# [55] tools_3.3.2                bit64_0.9-7                DESeq2_1.14.1             
# [58] survival_2.41-3            colorspace_1.3-2           cluster_2.0.6             
# [61] GenomicRanges_1.26.4       memoise_1.1.0              knitr_1.20             

