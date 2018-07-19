# R script for Gene Ontology analyses #


This script refers to the work performed by Severo et al, PNAS 2018 on single cell RNA Sequencing of Anopheles gambiae hemocytes. It includes the approach used to obtain the Gene Ontology (GO) analyses available in the above mentioned publication. 


## Getting started 

Prior to starting this script, you will need to obtain the mapped read counts for the scRNA-seq. The raw reads have been deposited under the accession number: "PRJEB23372" and instructions on how to perform the mapping are listed in the Materials and Methods section of the publication. 


### Installing packages

To run the R script described below, you will need to download and install the following packages either from CRAN or Bioconductor.:

install.packages(“plyr”)

source("https://bioconductor.org/biocLite.R")

biocLite(c("DESeq2", “topGO “, “genefilter”, "org.Ag.eg.db"))


## Running the script part I – preparing the data

A total of 26 cells and two pools of 30 cells each were sequenced in the work described at the referenced publication. In the first part of the script, you will prepare the data to be used for the Gene ontology (GO) analyses and normalize it based on the size factor of the ERCC spike-ins added for normalization and technical noise estimation purposes prior to cDNA synthesis using DESeq2 [1].


### to load mapped read counts, prepare the data and subset cells 

counts <- read.table("~/Desktop/GitHub/20170126-count.header.star.tab", header=T,  row.names=1)

columns <- c("A1", "A2", "A4", "A5", "A6", "A7", "B8", "C1", "C3", "C7", "C8", "D1", "D3", "D4", "D5", "D9", "E4", "F1", "F3", "F5", "F6", "F7", "F9", "F10", "G1", "G3", "G4", "H4")
all_counts <- subset(counts, select=columns)

features <- c("N_unmapped", "N_multimapping", "N_noFeature", "N_ambiguous")

all_counts_genes <- all_counts[-which(rownames(all_counts) %in% features), ]

all_counts_features <- all_counts[which(rownames(all_counts) %in% features), ]


### discard outlier cell (as discussed in the publication and indicated in Figure S1) and subset single cells from pools

sc_counts <- subset(all_counts_genes, select=-c(A1, A2, D1))


### to subset endogenous counts and ERCC spike-ins counts

geneTypes <- factor( c( AG="Ag", ER="ERCC" )[
  substr( rownames(sc_counts), 1, 2 ) ] )
countsAg <- sc_counts[ which( geneTypes=="Ag" ), ]
countsERCC <- sc_counts[ which( geneTypes=="ERCC" ), ]


### to divide by the size factors to get normalized counts (by sfERCC)

library("DESeq2")

sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
nCountsERCC <- t( t(countsERCC) / sfERCC )
nCountsAg <- t( t(countsAg) / sfERCC)


## Running the script part II - GO analyses using topGO and core transcriptome

The data is now normalized and ready for the second part of the script, where first you will filter the gene counts present in 90% of the cells to explore the existence of a “core” transcriptome and the GO terms that characterize these genes. You will use the topGO package [2] and the org.Ag.eg.db annotation package [3] and you will run both “classic” and “elim” methods of GO analyses for comparison purposes. The background gene set will be all genes present in the transcriptome.


### load the packages needed

library(topGO)
library(plyr)
library(genefilter)
source("https://bioconductor.org/biocLite.R")
biocLite("org.Ag.eg.db")


### Filter the genes present in 90% of the cells (Table S4)

sc_nCountsAg1 <- nCountsAg[rowSums(nCountsAg) >=1,]
filter_sc_cells <- apply(sc_nCountsAg1, 1, function(x) length(x[x>=1])>=22)
all_sc_cells <- sc_nCountsAg1 [filter_sc_cells,]


### to run GO analyses using topGO and 90% of the cells

geneNames <- rownames(nCountsAg)
geneList <- factor(as.integer(geneNames %in% rownames(all_sc_cells)))
names(geneList) <- geneNames
str(geneList)

GOdata_all_sc_cells_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Ag.eg.db", ID ="ensembl")

resultTopGO.elim_all_sc_cells_BP <- runTest(GOdata_all_sc_cells_BP, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic_all_sc_cells_BP <- runTest(GOdata_all_sc_cells_BP, algorithm = "classic", statistic = "Fisher" )

topGOres.elim_all_sc_cells_BP <- GenTable(GOdata_all_sc_cells_BP, classic = resultTopGO.elim_all_sc_cells_BP , orderBy = "classic", topNodes=200)
topGOres.clas_all_sc_cells_BP <- GenTable(GOdata_all_sc_cells_BP, classic = resultTopGO.classic_all_sc_cells_BP , orderBy = "classic", topNodes=200)



## Running the script part III - GO analyses using topGO and PPO high and PPO low populations

In this part, you will divide the data based on the identified subpopulations described in the manuscript. These populations were obtained using the technical noise estimation method described by Brennecke et al [4] followed by PCA analyses of the resulting variable genes. The R script used for that is available in the reference article and any modifications were listed in the Materials and Methods section of the Severo et al publication. The populations were compared based of differentially expressed (DE) genes obtained by using DESeq2.


### to set the groups that will be compared

samplenames2 <-c("Group 1", "Group 2",  "Group 1", "Group 1",	"Group 1",	"Group 2",	"Group 1", "Group 1",	"Group 2",	"Group 1", "Group 1", "Group 1", "Group 2",	"Group 1", "Group 1", "Group 1", "Group 1",	"Group 1", "Group 1", "Group 2", "Group 1", "Group 2", "Group 1")


### to filter out samples that are outliers on the PCA (Figure 2 of the publication)

countsAg2 <- subset(countsAg, select=-c(A4, B8))


### to estimate size factors without the outlier samples

countsERCC2 <- subset(countsERCC, select=-c(A4, B8))
sfERCC2 <- estimateSizeFactorsForMatrix( countsERCC2 )


### to build the object and run the DESeq2

conds = data.frame(samplenames2)
colnames(conds)="condition"
cds <- DESeqDataSetFromMatrix(countData = countsAg2, colData = conds, design = ~ condition)
sizeFactors(cds) <- sfERCC2
dds <- DESeq(cds, minReplicatesForReplace=Inf, fitType="local")
res=results( dds, contrast=c("condition","Group 2", "Group 1"), cooksCutoff=FALSE)
res<-res[order(res$padj),]


### to build a dataframe with the DE genes and inspect expression of low padj ones

res.df <- as.data.frame(res)
Deseq.genes <- res.df[res.df$padj <= 0.1,]
Deseq.genes <-na.omit(Deseq.genes)


### to make a data frame of the logFC results

Deseq.genes.pos <- Deseq.genes[Deseq.genes$log2FoldChange >= 0,]
Deseq.genes.neg <- Deseq.genes[Deseq.genes$log2FoldChange <= 0,]


### to compare the PPO populations based on GO terms (Table S6)

Below you will use the DE genes characterizing each population as input for topGO. The background gene set will be all genes present in the transcriptome.


#### to run topGO - Deseq pos

geneNames <- rownames(nCountsAg)
geneList <- factor(as.integer(geneNames %in% rownames(Deseq.genes.pos)))
names(geneList) <- geneNames
str(geneList)

GOdata_Deseq.genes.pos_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Ag.eg.db", ID ="ensembl")

resultTopGO.elim_Deseq.genes.pos_BP <- runTest(GOdata_Deseq.genes.pos_BP, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic_Deseq.genes.pos_BP <- runTest(GOdata_Deseq.genes.pos_BP, algorithm = "classic", statistic = "Fisher" )

topGOres.elim_Deseq.genes.pos_BP <- GenTable(GOdata_Deseq.genes.pos_BP, classic = resultTopGO.elim_Deseq.genes.pos_BP , orderBy = "classic", topNodes=200)
topGOres.clas_Deseq.genes.pos_BP <- GenTable(GOdata_Deseq.genes.pos_BP, classic = resultTopGO.classic_Deseq.genes.pos_BP , orderBy = "classic", topNodes=200)

#### to run topGO - Deseq neg

geneNames <- rownames(nCountsAg)
geneList <- factor(as.integer(geneNames %in% rownames(Deseq.genes.neg)))
names(geneList) <- geneNames
str(geneList)

GOdata_Deseq.genes.neg_BP <- new("topGOdata", ontology = "BP", allGenes = geneList, nodeSize=5, annot=annFUN.org, mapping="org.Ag.eg.db", ID ="ensembl")

resultTopGO.elim_Deseq.genes.neg_BP <- runTest(GOdata_Deseq.genes.neg_BP, algorithm = "elim", statistic = "Fisher" )
resultTopGO.classic_Deseq.genes.neg_BP <- runTest(GOdata_Deseq.genes.neg_BP, algorithm = "classic", statistic = "Fisher" )

topGOres.elim_Deseq.genes.neg_BP <- GenTable(GOdata_Deseq.genes.neg_BP, classic = resultTopGO.elim_Deseq.genes.neg_BP , orderBy = "classic", topNodes=200)
topGOres.clas_Deseq.genes.neg_BP <- GenTable(GOdata_Deseq.genes.neg_BP, classic = resultTopGO.classic_Deseq.genes.neg_BP , orderBy = "classic", topNodes=200)


## Author: 

Maiara Severo


## Publication: 

Severo MS, et al. Unbiased classification of mosquito blood cells by single-cell genomics and high-content imaging. Proceedings of the National Academy of Sciences of the United States of America. (in press)


## References:

1- Love MI, Huber W, & Anders S (2014) Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2. Genome Biol 15(12):550.

2- Alexa A, Rahnenfuhrer J, & Lengauer T (2006) Improved scoring of functional groups from gene expression data by decorrelating GO graph structure. Bioinformatics 22(13):1600-1607.

3 -Carlson M (2018). org.Ag.eg.db: Genome wide annotation for Anopheles.

4- Brennecke P, et al. (2013) Accounting for technical noise in single-cell RNA-seq experiments. Nat Methods 10(11):1093-1095.

