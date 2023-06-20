## 1. loading the required libraries

library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(EDASeq)
library(graphite)
library(houseOfClipUtility)
library(MOSClip)
library(curatedTCGAData)
library(TCGAutils)
library(tidyverse)
library(dplyr)
library(Hmisc)
library(readxl)
library(data.table)
library(gridExtra)
library(ggVennDiagram)
library(ggplot2)
library(reshape2)
library(igraph)
library(openxlsx)
library(penalized)
library(gridExtra)
library(cowplot)
library(survival)
library(survminer)

## 2. Pathway detailed analysis

### 2.1 Finding the pathways shared among different types of cancers

#Plotting the number of shared pathways 
BRCA <- rownames(pathwaySummaryHuge_BRCA)
CESC <- rownames(pathwaySummaryHuge_CESC)
OV <- rownames(pathwaySummaryHuge_OV)
UCEC <- rownames(pathwaySummaryHuge_UCEC)
UCS <- rownames(pathwaySummaryHuge_UCS)
x <- list(BRCA=BRCA, OV=OV, UCEC=UCEC, UCS=UCS, CESC=CESC)
ggVennDiagram(x, label = "count", label_alpha  = 0, )  + scale_fill_distiller(palette = "Pastel2")

#Plotting the percentage of pathways in common/unique
BRCA <- rownames(pathwaySummaryHuge_BRCA)
CESC <- rownames(pathwaySummaryHuge_CESC)
OV <- rownames(pathwaySummaryHuge_OV)
UCEC <- rownames(pathwaySummaryHuge_UCEC)
UCS <- rownames(pathwaySummaryHuge_UCS)
x <- list(BRCA=BRCA, OV=OV, UCEC=UCEC, UCS=UCS, CESC=CESC)
ggVennDiagram(x, label = "percent", label_alpha  = 0, label_percent_digit = 1)  + scale_fill_distiller(palette = "Set3")

### 2.2 Checking the names of the pathways that were in common

#OV_UCEC_CESC = 4 pathways
OV_UCEC <- intersect(rownames(pathwaySummaryHuge_OV), rownames(pathwaySummaryHuge_UCEC))
OV_UCEC_CESC <- intersect(OV_UCEC, rownames(pathwaySummaryHuge_CESC))

#BRCA_UCEC_OV = 3 pathways
BRCA_UCEC <- intersect(rownames(pathwaySummaryHuge_BRCA), rownames(pathwaySummaryHuge_UCEC))
BRCA_UCEC_OV <- intersect(BRCA_UCEC, rownames(pathwaySummaryHuge_OV))

#CESC_OV_BRCA = 1 pathway
CESC_OV <- intersect(rownames(pathwaySummaryHuge_CESC), rownames(pathwaySummaryHuge_OV))
CESC_OV_BRCA <- intersect(CESC_OV, rownames(pathwaySummaryHuge_BRCA))

#BRCA_CESC_UCEC = 24 pathways
BRCA_CESC <- intersect(rownames(pathwaySummaryHuge_BRCA), rownames(pathwaySummaryHuge_CESC))
BRCA_CESC_UCEC <- intersect(BRCA_CESC, rownames(pathwaySummaryHuge_UCEC))

# Plotting the whole groups of cancers with their pathways
# Create a data frame
data <- data.frame(
  variable = rep(c("OV_UCEC_CESC", "BRCA_UCEC_OV", "CESC_OV_BRCA", "BRCA_CESC_UCEC"), times = c(length(OV_UCEC_CESC), length(BRCA_UCEC_OV), length(CESC_OV_BRCA), length(BRCA_CESC_UCEC))),
  element = c(OV_UCEC_CESC, BRCA_UCEC_OV, CESC_OV_BRCA, BRCA_CESC_UCEC)
)

# Convert element column to factor and specify levels
data$element <- factor(data$element, levels = unique(data$element))

# Manually specify a color palette with distinct colors
colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
            "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494",
            "#B3B3B3", "#8DD3C7", "#FB8072", "#80B1D3", "#FDB462", "#BEBADA", "#FCCDE5", "#D9D9D9",
            "#BC80BD", "#CCEBC5", "#FFED6F", "#7FC97F", "#BEAED4", "#FDC086", "#FFFF99", "#386CB0")

# Create a stacked bar chart
p <- ggplot(data, aes(x = variable, fill = element)) +
  geom_bar() +
  scale_fill_manual(values = colors) +
  theme(legend.position = "bottom") +
  labs(fill = "Pathways")
ggsave("plot.png", plot = p, width = 23, height = 10)

## 3. Filtering the edges by each module

### 3.1 BRCA

df_BRCA <- data.frame()
n <- 0
for (pathway in unique(moduleSummary_BRCA$pathway)){
    n <- n + 1
    print(pathway)
    graphnel <- pathwayGraph(reactome[[pathway]], which = "proteins")
    igraph <- igraph::igraph.from.graphNEL(graphnel)
    genes <- list()
    for (module in rownames(moduleSummary_BRCA[moduleSummary_BRCA$pathway == pathway,])){
      str <- module
      a <- strsplit(str, "[.]")
      m <- a[[1]][2]
      print(module)
      genes <- c(genes, multiOmics_module_BRCA[[pathway]]@modules[[as.integer(m)]])
    }
    genes <- as.character(genes)
    genes <- genes[!duplicated(genes)]
    subgraph <- subgraph(igraph, genes)
    df_BRCA[n, "pathway"] <- pathway
    df_BRCA[n, "ratio"] <- (gsize(subgraph) / length(V(subgraph)))
    print(gsize(subgraph))
    print(length(V(subgraph)))
}

plot(df_BRCA$ratio, type = "b", xlab = "Graphs", xaxt = "n", ylab = "Average Number of Edges for Each Node") 

boxplot(df_BRCA$ratio, ylab = "Average Number of Edges for Each Node")

### 3.2 CESC

df_CESC <- data.frame()
n <- 0
for (pathway in unique(moduleSummary_CESC$pathway)){
    n <- n + 1
    print(pathway)
    graphnel <- pathwayGraph(reactome[[pathway]], which = "proteins")
    igraph <- igraph::igraph.from.graphNEL(graphnel)
    genes <- list()
    for (module in rownames(moduleSummary_CESC[moduleSummary_CESC$pathway == pathway,])){
      str <- module
      a <- strsplit(str, "[.]")
      m <- a[[1]][2]
      print(module)
      genes <- c(genes, multiOmics_module_CESC[[pathway]]@modules[[as.integer(m)]])
    }
    genes <- as.character(genes)
    genes <- genes[!duplicated(genes)]
    subgraph <- subgraph(igraph, genes)
    df_CESC[n, "pathway"] <- pathway
    df_CESC[n, "ratio"] <- (gsize(subgraph) / length(V(subgraph)))
    print(gsize(subgraph))
    print(length(V(subgraph)))
}

plot(df_CESC$ratio, type = "b", xlab = "Graphs", xaxt = "n", ylab = "Average Number of Edges for Each Node") 

boxplot(df_CESC$ratio, ylab = "Average Number of Edges for Each Node")

### 3.3 OV

df_OV <- data.frame()
n <- 0
for (pathway in unique(moduleSummary_OV$pathway)){
    n <- n + 1
    print(pathway)
    graphnel <- pathwayGraph(reactome[[pathway]], which = "proteins")
    igraph <- igraph::igraph.from.graphNEL(graphnel)
    genes <- list()
    for (module in rownames(moduleSummary_OV[moduleSummary_OV$pathway == pathway,])){
      str <- module
      a <- strsplit(str, "[.]")
      m <- a[[1]][2]
      print(module)
      genes <- c(genes, multiOmics_module_OV[[pathway]]@modules[[as.integer(m)]])
    }
    genes <- as.character(genes)
    genes <- genes[!duplicated(genes)]
    subgraph <- subgraph(igraph, genes)
    df_OV[n, "pathway"] <- pathway
    df_OV[n, "ratio"] <- (gsize(subgraph) / length(V(subgraph)))
    print(gsize(subgraph))
    print(length(V(subgraph)))
}

plot(df_OV$ratio, type = "b", xlab = "Graphs", xaxt = "n", ylab = "Average Number of Edges for Each Node") 

boxplot(df_OV$ratio, ylab = "Average Number of Edges for Each Node")

### 3.4 UCEC

df_UCEC <- data.frame()
n <- 0
for (pathway in unique(moduleSummary_UCEC$pathway)){
    n <- n + 1
    print(pathway)
    graphnel <- pathwayGraph(reactome[[pathway]], which = "proteins")
    igraph <- igraph::igraph.from.graphNEL(graphnel)
    genes <- list()
    for (module in rownames(moduleSummary_UCEC[moduleSummary_UCEC$pathway == pathway,])){
      str <- module
      a <- strsplit(str, "[.]")
      m <- a[[1]][2]
      print(module)
      genes <- c(genes, multiOmics_module_UCEC[[pathway]]@modules[[as.integer(m)]])
    }
    genes <- as.character(genes)
    genes <- genes[!duplicated(genes)]
    subgraph <- subgraph(igraph, genes)
    df_UCEC[n, "pathway"] <- pathway
    df_UCEC[n, "ratio"] <- (gsize(subgraph) / length(V(subgraph)))
    print(gsize(subgraph))
    print(length(V(subgraph)))
}

plot(df_UCEC$ratio, type = "b", xlab = "Graphs", xaxt = "n", ylab = "Average Number of Edges for Each Node") 

boxplot(df_UCEC$ratio, ylab = "Average Number of Edges for Each Node")

### 3.5 UCS

df_UCS <- data.frame()
n <- 0
for (pathway in unique(moduleSummary_UCS$pathway)){
    n <- n + 1
    print(pathway)
    graphnel <- pathwayGraph(reactome[[pathway]], which = "proteins")
    igraph <- igraph::igraph.from.graphNEL(graphnel)
    genes <- list()
    for (module in rownames(moduleSummary_UCS[moduleSummary_UCS$pathway == pathway,])){
      str <- module
      a <- strsplit(str, "[.]")
      m <- a[[1]][2]
      print(module)
      genes <- c(genes, multiOmics_module_UCS[[pathway]]@modules[[as.integer(m)]])
    }
    genes <- as.character(genes)
    genes <- genes[!duplicated(genes)]
    subgraph <- subgraph(igraph, genes)
    df_UCS[n, "pathway"] <- pathway
    df_UCS[n, "ratio"] <- (gsize(subgraph) / length(V(subgraph)))
    print(gsize(subgraph))
    print(length(V(subgraph)))
}

plot(df_UCS$ratio, type = "b", xlab = "Graphs", xaxt = "n", ylab = "Average Number of Edges for Each Node") 

boxplot(df_UCS$ratio, ylab = "Average Number of Edges for Each Node")

## 4. Preparing the matrix for the Cytoscape

cytoscape <- function(moduleSummary, cancer_name, multiOmics){
  cyto_df <- data.frame()
  for (module in rownames(moduleSummary)){
    print(module)
    genes <- list()
    omics <- list()
    str <- module
    a <- strsplit(str, "[.]")
    pathway <- a[[1]][1]
    m <- a[[1]][2]
    graphnel <- pathwayGraph(reactome[[pathway]], which = "proteins")
    igraph <- igraph::igraph.from.graphNEL(graphnel)
    genes <- multiOmics[[pathway]]@modules[[as.integer(m)]]
    genes <- as.character(genes)
    genes <- genes[!duplicated(genes)]
    subgraph <- subgraph(igraph, genes)
    #We can use this condition to filter the modules with too many edges
    if ((gsize(subgraph) / length(V(subgraph))) <= 40 & (gsize(subgraph) / length(V(subgraph))) > 0){
      cyto_df_new <- data.frame()
      cyto_df_new <- data.frame(igraph::as_edgelist(subgraph))
      cyto_df_new$cancer <- cancer_name
      cyto_df_new$pathway <- pathway
      cyto_df_new$module <- m
      for (col in 4:ncol(moduleSummary)) {
        if (moduleSummary[module, col] <= 0.05 & !is.na(moduleSummary[module, col])) {
          omics <- c(omics, substr(colnames(moduleSummary)[col], 1, 3))
        }
      omics <- unique(unlist(omics))
      cyto_df_new$omics <- paste(omics, collapse = "-")
      }
      cyto_df <- rbind(cyto_df, cyto_df_new)
    }
  }
  return(cyto_df)
}

### 4.1 BRCA

cyto_BRCA <- cytoscape(moduleSummary = moduleSummary_BRCA, cancer_name = "BRCA", multiOmics = multiOmics_module_BRCA)

### 4.2 CESC

cyto_CESC <- cytoscape(moduleSummary = moduleSummary_CESC, cancer_name = "CESC", multiOmics = multiOmics_module_CESC)

### 4.3 OV

cyto_OV <- cytoscape(moduleSummary = moduleSummary_OV, cancer_name = "OV", multiOmics = multiOmics_module_OV)

### 4.4 UCEC

cyto_UCEC <- cytoscape(moduleSummary = moduleSummary_UCEC, cancer_name = "UCEC", multiOmics = multiOmics_module_UCEC)

### 4.5 UCS

cyto_UCS <- cytoscape(moduleSummary = moduleSummary_UCS, cancer_name = "UCS", multiOmics = multiOmics_module_UCS)

## 5. Conveting the gene ids

### 5.1 BRCA

#This step must be repeated for both X1 and X2 columns
for (gene in unique(cyto_BRCA$X2)){
  a <- strsplit(gene, ":")
  cyto_BRCA[cyto_BRCA$X2 == gene, "X2"]<- a[[1]][2]
}
symbol_X2 <- AnnotationDbi::select(org.Hs.eg.db, keys= cyto_BRCA$X2, keytype = "ENTREZID" ,columns = "SYMBOL")
cyto_BRCA$X2 <- symbol_X2$SYMBOL

### 5.2 CESC

for (gene in unique(cyto_CESC$X1)){
  a <- strsplit(gene, ":")
  cyto_CESC[cyto_CESC$X1 == gene, "X1"]<- a[[1]][2]
}
symbol_X1 <- AnnotationDbi::select(org.Hs.eg.db, keys= cyto_CESC$X1, keytype = "ENTREZID" ,columns = "SYMBOL")
cyto_CESC$X1 <- symbol_X1$SYMBOL

### 5.3 OV

for (gene in unique(cyto_OV$X2)){
  a <- strsplit(gene, ":")
  cyto_OV[cyto_OV$X2 == gene, "X2"]<- a[[1]][2]
}
symbol_X2 <- AnnotationDbi::select(org.Hs.eg.db, keys= cyto_OV$X2, keytype = "ENTREZID" ,columns = "SYMBOL")
cyto_OV$X2 <- symbol_X2$SYMBOL

### 5.4 UCEC

for (gene in unique(cyto_UCEC$X1)){
  a <- strsplit(gene, ":")
  cyto_UCEC[cyto_UCEC$X1 == gene, "X1"]<- a[[1]][2]
}
symbol_X1 <- AnnotationDbi::select(org.Hs.eg.db, keys= cyto_UCEC$X1, keytype = "ENTREZID" ,columns = "SYMBOL")
cyto_UCEC$X1 <- symbol_X1$SYMBOL

### 5.5 UCS

for (gene in unique(cyto_UCS$X2)){
  a <- strsplit(gene, ":")
  cyto_UCS[cyto_UCS$X2 == gene, "X2"]<- a[[1]][2]
}
symbol_X2 <- AnnotationDbi::select(org.Hs.eg.db, keys= cyto_UCS$X2, keytype = "ENTREZID" ,columns = "SYMBOL")
cyto_UCS$X2 <- symbol_X2$SYMBOL

## 6. Creating the plot for all genes together

### 6.1 BRCA

for (i in c(1:12040)){
  cyto_BRCA[(i + 12040), "X1"] <- cyto_BRCA$X2[i]
  cyto_BRCA[(i + 12040), "X2"] <- cyto_BRCA$X1[i]
  cyto_BRCA[(i + 12040), "cancer"] <- "BRCA"
  cyto_BRCA[(i + 12040), "pathway"] <- cyto_BRCA$pathway[i]
  cyto_BRCA[(i + 12040), "module"] <- cyto_BRCA$module[i]
  cyto_BRCA[(i + 12040), "omics"] <- cyto_BRCA$omics[i]
}

save(cyto_BRCA, file = "./cyto_BRCA.RData")
write.xlsx(cyto_BRCA, "./cyto_BRCA.xls", sep = "\t")

### 6.2 CESC

for (i in c(1:471)){
  cyto_CESC[(i + 471), "X1"] <- cyto_CESC$X2[i]
  cyto_CESC[(i + 471), "X2"] <- cyto_CESC$X1[i]
  cyto_CESC[(i + 471), "cancer"] <- "CESC"
  cyto_CESC[(i + 471), "pathway"] <- cyto_CESC$pathway[i]
  cyto_CESC[(i + 471), "module"] <- cyto_CESC$module[i]
  cyto_CESC[(i + 471), "omics"] <- cyto_CESC$omics[i]
}

save(cyto_CESC, file = "./cyto_CESC.RData")
write.xlsx(cyto_CESC, "./cyto_CESC.xls", sep = "\t")

### 6.3 OV

for (i in c(1:92)){
  cyto_OV[(i + 92), "X1"] <- cyto_OV$X2[i]
  cyto_OV[(i + 92), "X2"] <- cyto_OV$X1[i]
  cyto_OV[(i + 92), "cancer"] <- "OV"
  cyto_OV[(i + 92), "pathway"] <- cyto_OV$pathway[i]
  cyto_OV[(i + 92), "module"] <- cyto_OV$module[i]
  cyto_OV[(i + 92), "omics"] <- cyto_OV$omics[i]
}

save(cyto_OV, file = "./cyto_OV.RData")
write.xlsx(cyto_OV, "./cyto_OV.xls", sep = "\t")

### 6.4 UCEC

for (i in c(1:7372)){
  cyto_UCEC[(i + 7372), "X1"] <- cyto_UCEC$X2[i]
  cyto_UCEC[(i + 7372), "X2"] <- cyto_UCEC$X1[i]
  cyto_UCEC[(i + 7372), "cancer"] <- "UCEC"
  cyto_UCEC[(i + 7372), "pathway"] <- cyto_UCEC$pathway[i]
  cyto_UCEC[(i + 7372), "module"] <- cyto_UCEC$module[i]
  cyto_UCEC[(i + 7372), "omics"] <- cyto_UCEC$omics[i]
}

save(cyto_UCEC, file = "./cyto_UCEC.RData")
write.xlsx(cyto_UCEC, "./cyto_UCEC.xls", sep = "\t")

### 6.5 UCS

for (i in c(1:185)){
  cyto_UCS[(i + 185), "X1"] <- cyto_UCS$X2[i]
  cyto_UCS[(i + 185), "X2"] <- cyto_UCS$X1[i]
  cyto_UCS[(i + 185), "cancer"] <- "UCS"
  cyto_UCS[(i + 185), "pathway"] <- cyto_UCS$pathway[i]
  cyto_UCS[(i + 185), "module"] <- cyto_UCS$module[i]
  cyto_UCS[(i + 185), "omics"] <- cyto_UCS$omics[i]
}

save(cyto_UCS, file = "./cyto_UCS.RData")
write.xlsx(cyto_UCS, "./cyto_UCS.xls", sep = "\t")

## 7. Extracted all of the significant genes from each cancer

gene_list <- function(moduleSummary, multiOmics_module){
  genes <- list()
  for (module in rownames(moduleSummary)){
    a <- strsplit(module, "[.]")
    pathway <- a[[1]][1]
    m <- a[[1]][2]
    genes <- c(genes, multiOmics_module[[pathway]]@modules[[as.integer(m)]])
  }
  genes <- as.character(genes)
  genes <- genes[!duplicated(genes)]
  return(genes)
}

### 7.1 BRCA

BRCA_genes <- gene_list(moduleSummary_BRCA, multiOmics_module_BRCA)

### 7.2 CESC

CESC_genes <- gene_list(moduleSummary_CESC, multiOmics_module_CESC)

### 7.3 OV

OV_genes <- gene_list(moduleSummary_OV, multiOmics_module_OV)

### 7.4 UCEC

UCEC_genes <- gene_list(moduleSummary_UCEC, multiOmics_module_UCEC)

### 7.5 UCS

UCS_genes <- gene_list(moduleSummary_UCS, multiOmics_module_UCS)

## 8. Making a dataframe for each omic

make_matrix <- function(omic_df, genes, type){
  matrix <- data.frame()
  a <- omic_df[genes,]
  a <- na.omit(a)
  for (i in c(1:(dim(a)[1]))){
    rownames(a)[i] <- paste0(rownames(a)[i], ".", type)
  }
  return(a)
}

### 8.1 BRCA

#### 8.1.1 Expression

matrix_BRCA_exp <- make_matrix(BRCA_EXP, BRCA_genes, "exp")
matrix_BRCA_exp[] <- lapply(matrix_BRCA_exp, as.numeric)

#### 8.1.2 CNV

matrix_BRCA_cnv <- make_matrix(BRCA_CNV, BRCA_genes, "cnv")
matrix_BRCA_cnv[] <- lapply(matrix_BRCA_cnv, as.numeric)

#### 8.1.3 Methylation

matrix_BRCA_met <- make_matrix(BRCA_Methylation, BRCA_genes, "met")
matrix_BRCA_met[] <- lapply(matrix_BRCA_met, as.numeric)

#### 8.1.4 Mutaion

matrix_BRCA_mut <- make_matrix(BRCA_Mutation, BRCA_genes, "mut")
matrix_BRCA_mut[] <- lapply(matrix_BRCA_mut, as.numeric)

### 8.2 CESC

#### 8.2.1 Expression

matrix_CESC_exp <- make_matrix(CESC_EXP, CESC_genes, "exp")
matrix_CESC_exp[] <- lapply(matrix_CESC_exp, as.numeric)

#### 8.2.2 CNV

matrix_CESC_cnv <- make_matrix(CESC_CNV, CESC_genes, "cnv")
matrix_CESC_cnv[] <- lapply(matrix_CESC_cnv, as.numeric)

#### 8.2.3 Methylation

matrix_CESC_met <- make_matrix(CESC_Methylation, CESC_genes, "met")
matrix_CESC_met[] <- lapply(matrix_CESC_met, as.numeric)

#### 8.2.4 Mutaion

matrix_CESC_mut <- make_matrix(CESC_Mutation, CESC_genes, "mut")
matrix_CESC_mut[] <- lapply(matrix_CESC_mut, as.numeric)

### 8.3 OV

#### 8.3.1 Expression

matrix_OV_exp <- make_matrix(OV_EXP, OV_genes, "exp")
matrix_OV_exp[] <- lapply(matrix_OV_exp, as.numeric)

#### 8.3.2 CNV

matrix_OV_cnv <- make_matrix(OV_CNV, OV_genes, "cnv")
matrix_OV_cnv[] <- lapply(matrix_OV_cnv, as.numeric)

#### 8.3.3 Methylation

matrix_OV_met <- make_matrix(OV_Methylation, OV_genes, "met")
matrix_OV_met[] <- lapply(matrix_OV_met, as.numeric)

#### 8.3.4 Mutaion

matrix_OV_mut <- make_matrix(OV_Mutation, OV_genes, "mut")
matrix_OV_mut[] <- lapply(matrix_OV_mut, as.numeric)

### 8.4 UCEC

#### 8.4.1 Expression

matrix_UCEC_exp <- make_matrix(UCEC_EXP, UCEC_genes, "exp")
matrix_UCEC_exp[] <- lapply(matrix_UCEC_exp, as.numeric)

#### 8.4.2 CNV

matrix_UCEC_cnv <- make_matrix(UCEC_CNV, UCEC_genes, "cnv")
matrix_UCEC_cnv[] <- lapply(matrix_UCEC_cnv, as.numeric)

#### 8.4.3 Methylation

matrix_UCEC_met <- make_matrix(UCEC_Methylation, UCEC_genes, "met")
matrix_UCEC_met[] <- lapply(matrix_UCEC_met, as.numeric)

#### 8.4.4 Mutaion

matrix_UCEC_mut <- make_matrix(UCEC_Mutation, UCEC_genes, "mut")
matrix_UCEC_mut[] <- lapply(matrix_UCEC_mut, as.numeric)

### 8.5 UCS

#### 8.5.1 Expression

matrix_UCS_exp <- make_matrix(UCS_EXP, UCS_genes, "exp")
matrix_UCS_exp[] <- lapply(matrix_UCS_exp, as.numeric)

#### 8.5.2 CNV

matrix_UCS_cnv <- make_matrix(UCS_CNV, UCS_genes, "cnv")
matrix_UCS_cnv[] <- lapply(matrix_UCS_cnv, as.numeric)

#### 8.5.3 Methylation

matrix_UCS_met <- make_matrix(UCS_Methylation, UCS_genes, "met")
matrix_UCS_met[] <- lapply(matrix_UCS_met, as.numeric)

#### 8.5.4 Mutaion

matrix_UCS_mut <- make_matrix(UCS_Mutation, UCS_genes, "mut")
matrix_UCS_mut[] <- lapply(matrix_UCS_mut, as.numeric)

## 9. Writing a function to calculate the 0 and 1 for each gene/patient

#A function for expression and methylation
exp_met_matrix <- function(survival, matrix, l){
  penalized <- penalized::penalized(Surv(survival$days, survival$status), penalized = as.data.frame(t(matrix)), unpenalized = ~0, lambda1 = l, model = "cox")
  HR <- as.data.frame(coefficients(penalized))
  matrix <- matrix[rownames(matrix) %in% rownames(HR), ]
  for (gene in rownames(matrix)){
    median <- median(as.numeric(matrix[gene,]))
    values <- as.numeric(matrix[gene,])
    max <- (max(values) + 1)
    if (HR[gene,] < 0){
      values[values >= median] <- max
      values[values < median] <- abs(as.numeric(HR[gene,]))
      values[values == max] <- 0
    }  
    else if(HR[gene,] > 0){
      values[values > median] <- max
      values[values <= median] <-  0
      values[values == max] <- abs(as.numeric(HR[gene,]))
    }
    matrix[gene,] <- values
  }
  return(matrix)
}

#A function for CNV
cnv_matrix <- function(survival, matrix, l){
  penalized <- penalized::penalized(Surv(survival$days, survival$status), penalized = as.data.frame(t(matrix)), unpenalized = ~0, lambda1 = l, model = "cox")
  HR <- as.data.frame(coefficients(penalized))
  matrix <- matrix[rownames(matrix) %in% rownames(HR), ]
  for (gene in rownames(matrix)){
    values <- as.numeric(matrix[gene,])
    if (HR[gene,] < 0){
      values[values >= 0] <- 0
      values[values < 0] <- abs(as.numeric(HR[gene,]))
    }
    else if(HR[gene,] > 0){
      values[values > 0] <- abs(as.numeric(HR[gene,]))
      values[values <= 0] <- 0
    }
    matrix[gene,] <- values
  }
  return(matrix)
}

#A function for mutation
mutation_matrix <- function(survival, matrix, l){
  penalized <- penalized::penalized(Surv(survival$days, survival$status), penalized = as.data.frame(t(matrix)), unpenalized = ~0, lambda1 = l, model = "cox")
  HR <- as.data.frame(coefficients(penalized))
  matrix <- matrix[rownames(matrix) %in% rownames(HR), ]
  for (gene in rownames(matrix)){
    values <- as.numeric(matrix[gene,])
    if (HR[gene,] < 0){
      values[values == 1] <- 0
    }
    else if(HR[gene,] > 0){
      values[values == 1] <- abs(as.numeric(HR[gene,]))
    }
    matrix[gene,] <- values
  }
  return(matrix)
}

### 9.1 BRCA

#### 9.1.1 Expression

BRCA_exp_final <- exp_met_matrix(survival_BRCA, matrix_BRCA_exp, 3.4)

#### 9.1.2 Methylation

BRCA_met_final <- exp_met_matrix(survival_BRCA, matrix_BRCA_met, 3.4)

#### 9.1.3 CNV

BRCA_cnv_final <- cnv_matrix(survival_BRCA, matrix_BRCA_cnv, 3.4)

#### 9.1.4 Mutation

BRCA_mut_final <- mutation_matrix(survival_BRCA, matrix_BRCA_mut, 3.4)

matrix_BRCA <- do.call("rbind", list(BRCA_exp_final, BRCA_cnv_final, BRCA_met_final, BRCA_mut_final))
save(matrix_BRCA, file = "./matrix_BRCA.RData")

### 9.2 CESC

#### 9.2.1 Expression

CESC_exp_final <- exp_met_matrix(survival_CESC, matrix_CESC_exp, 1.4)

#### 9.2.2 Methylation

CESC_met_final <- exp_met_matrix(survival_CESC, matrix_CESC_met, 1.4)

#### 9.2.3 CNV

CESC_cnv_final <- cnv_matrix(survival_CESC, matrix_CESC_cnv, 1.4)

#### 9.2.4 Mutataion

CESC_mut_final <- mutation_matrix(survival_CESC, matrix_CESC_mut, 1.4)

matrix_CESC <- do.call("rbind", list(CESC_exp_final, CESC_cnv_final, CESC_met_final, CESC_mut_final))
save(matrix_CESC, file = "./matrix_CESC.RData")

### 9.3 OV

#### 9.3.1 Expression

OV_exp_final <- exp_met_matrix(survival_OV, matrix_OV_exp, 3.3)

#### 9.3.2 Methylation

OV_met_final <- exp_met_matrix(survival_OV, matrix_OV_met, 3.3)

#### 9.3.3 CNV

OV_cnv_final <- cnv_matrix(survival_OV, matrix_OV_cnv, 3.3)

#### 9.3.4 Mutataion

OV_mut_final <- mutation_matrix(survival_OV, matrix_OV_mut, 3.3)

matrix_OV <- do.call("rbind", list(OV_exp_final, OV_cnv_final, OV_met_final, OV_mut_final))
save(matrix_OV, file = "./matrix_OV.RData")

### 9.4 UCEC

#### 9.4.1 Expression

UCEC_exp_final <- exp_met_matrix(survival_UCEC, matrix_UCEC_exp, 3)

#### 9.4.2 Methylation

UCEC_met_final <- exp_met_matrix(survival_UCEC, matrix_UCEC_met, 3)

#### 9.4.3 CNV

UCEC_cnv_final <- cnv_matrix(survival_UCEC, matrix_UCEC_cnv, 3)

#### 9.4.4 Mutataion

UCEC_mut_final <- mutation_matrix(survival_UCEC, matrix_UCEC_mut, 3)

matrix_UCEC <- do.call("rbind", list(UCEC_exp_final, UCEC_cnv_final, UCEC_met_final, UCEC_mut_final))
save(matrix_UCEC, file = "./matrix_UCEC.RData")

### 9.5 UCS

#### 9.5.1 Expression

UCS_exp_final <- exp_met_matrix(survival_UCS, matrix_UCS_exp, 1.5)

#### 9.5.2 Methylation

UCS_met_final <- exp_met_matrix(survival_UCS, matrix_UCS_met, 1.5)

#### 9.5.3 CNV

UCS_cnv_final <- cnv_matrix(survival_UCS, matrix_UCS_cnv, 1.5)

#### 9.5.4 Mutataion

UCS_mut_final <- mutation_matrix(survival_UCS, matrix_UCS_mut, 1.5)

matrix_UCS <- do.call("rbind", list(UCS_exp_final, UCS_cnv_final, UCS_met_final, UCS_mut_final))
save(matrix_UCS, file = "./matrix_UCS.RData")

## 10. Creating data frames for survival heatmap and creating survival heatmaps 

heatmap_survival <- function(matrix){
  #changing the gene id to symbol
  genes <- list()
  for (gene in rownames(matrix)){
    b <- strsplit(gene, ":")
    c <- b[[1]][2]
    d <- strsplit(c, "[.]")
    genes <- c(genes, d[[1]][1])
  }
  genes <- as.character(genes)
  symbol <- AnnotationDbi::select(org.Hs.eg.db, keys= genes, keytype = "ENTREZID" ,columns = "SYMBOL")
  matrix$symbol <- symbol$SYMBOL
  matrix <- matrix %>% relocate(symbol)
  n <- 0
  for (gene in rownames(matrix)){
    n <- n+1
    d <- strsplit(gene, "[.]")
    matrix$symbol[n] <- paste0(matrix$symbol[n], ".", d[[1]][2])
  }
  rownames(matrix) <- matrix$symbol
  matrix$symbol <- NULL
  #making the total survival score
  patients_matrix <- colnames(matrix)
  for (patient in patients_matrix){
    matrix["SCORE", patient] <- 0
    matrix["SCORE", patient] <- sum(matrix[, patient])
  }
  matrix <- as.data.frame(t(matrix))
  matrix <- matrix[order(matrix$SCORE),]
  #rescaling the scores value
  q <- as.list(quantile(matrix$SCORE))
  matrix[matrix$SCORE <= q[[2]], "SCORE"] <- 0.2
  matrix[matrix$SCORE > q[[2]] & matrix$SCORE <= q[[4]], "SCORE"] <- 0.5
  matrix[matrix$SCORE > q[[4]] & matrix$SCORE <= q[[5]], "SCORE"] <- 0.8
  matrix <- matrix %>%
    rownames_to_column() %>%
    gather(colname, value, -rowname)
  colnames(matrix) <- c("patient", "gene.omic", "value")
  matrix$patient <- factor(matrix$patient, levels = unique(matrix$patient))
  matrix$gene.omic <- factor(matrix$gene.omic, levels = unique(matrix$gene.omic))
  matrix$value <- as.numeric(as.character(matrix$value))
  return(matrix)
  }

### 10.1 BRCA

# Increase all non-zero values by 1 to have better colors
matrix_BRCA[matrix_BRCA != 0] <- matrix_BRCA[matrix_BRCA != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_BRCA)
# Reorder the rows of the data frame based on their row sums
matrix_BRCA <- matrix_BRCA[order(row_sums, decreasing = FALSE), ]
heatmap_BRCA <- heatmap_survival(matrix_BRCA)

#Creating two heatmaps on top of each other
heatmap_BRCA$fill_color <- ifelse(heatmap_BRCA$gene.omic == "SCORE" & heatmap_BRCA$value == 0.2, "lightblue", ifelse(heatmap_BRCA$gene.omic == "SCORE" & heatmap_BRCA$value == 0.5, "orange", ifelse(heatmap_BRCA$gene.omic == "SCORE" & heatmap_BRCA$value == 0.8, "darkred", scales::col_numeric(c("white", "red"), domain = c(min(heatmap_BRCA$value), max(heatmap_BRCA$value)))(heatmap_BRCA$value))))

p1 <- ggplot(heatmap_BRCA[heatmap_BRCA$gene.omic == "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())

p2 <- ggplot(heatmap_BRCA[heatmap_BRCA$gene.omic != "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 2), axis.text.x = element_text(size = 0.5, angle = 90))

pdf("heatmap_BRCA.pdf")
plot_grid(p1,p2,nrow=2,align='v',axis='tb',rel_heights=c(1/20,19/20))
dev.off()

#legends for all types of cancers
heatmap_BRCA$fill_color <- ifelse(heatmap_BRCA$gene.omic == "SCORE" & heatmap_BRCA$value == 0.2, "lightblue", ifelse(heatmap_BRCA$gene.omic == "SCORE" & heatmap_BRCA$value == 0.5, "orange", ifelse(heatmap_BRCA$gene.omic == "SCORE" & heatmap_BRCA$value == 0.8, "darkred", scales::col_numeric(c("white", "red"), domain = c(min(heatmap_BRCA$value), max(heatmap_BRCA$value)))(heatmap_BRCA$value))))

p1 <- ggplot(heatmap_BRCA[heatmap_BRCA$gene.omic == "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())+ 
  scale_fill_manual(values = c("darkred", "orange", "lightblue"),
        name = "Risk of Death",
        labels = c("High", "Medium", "Low")) 

#legends
heatmap_BRCA$fill_color <- ifelse(heatmap_BRCA$gene.omic == "SCORE" & heatmap_BRCA$value == 0.2, "lightblue", ifelse(heatmap_BRCA$gene.omic == "SCORE" & heatmap_BRCA$value == 0.5, "orange", ifelse(heatmap_BRCA$gene.omic == "SCORE" & heatmap_BRCA$value == 0.8, "darkred", scales::col_numeric(c("white", "red"), domain = c(min(heatmap_BRCA$value), max(heatmap_BRCA$value)))(heatmap_BRCA$value))))

p2 <- ggplot(heatmap_BRCA[heatmap_BRCA$gene.omic != "SCORE", ], aes(x = patient, y = gene.omic, fill = value)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 2), axis.text.x = element_text(size = 0.5, angle = 90)) +
  scale_fill_gradient(low = "white", high = "red", name = "Hazard of Genes",
                      breaks  = c(min(heatmap_BRCA$value), max(heatmap_BRCA$value)/2, max(heatmap_BRCA$value) - 0.04),
                      labels = c("Low", "Intermediate", "High"))

### 10.2 CESC

# Increase all non-zero values by 1 to have better colors
matrix_CESC[matrix_CESC != 0] <- matrix_CESC[matrix_CESC != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_CESC)
# Reorder the rows of the data frame based on their row sums
matrix_CESC <- matrix_CESC[order(row_sums, decreasing = FALSE), ]
heatmap_CESC <- heatmap_survival(matrix_CESC)

#Creating two heatmaps on top of each other
heatmap_CESC$fill_color <- ifelse(heatmap_CESC$gene.omic == "SCORE" & heatmap_CESC$value == 0.2, "lightblue", ifelse(heatmap_CESC$gene.omic == "SCORE" & heatmap_CESC$value == 0.5, "orange", ifelse(heatmap_CESC$gene.omic == "SCORE" & heatmap_CESC$value == 0.8, "darkred", scales::col_numeric(c("white", "red"), domain = c(min(heatmap_CESC$value), max(heatmap_CESC$value)))(heatmap_CESC$value))))

p1 <- ggplot(heatmap_CESC[heatmap_CESC$gene.omic == "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())

p2 <- ggplot(heatmap_CESC[heatmap_CESC$gene.omic != "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 6), axis.text.x = element_text(size = 2, angle = 90))

pdf("heatmap_CESC.pdf")
plot_grid(p1,p2,nrow=2,align='v',axis='tb',rel_heights=c(1/20,19/20))
dev.off()

### 10.3 OV

# Increase all non-zero values by 1 to have better colors
matrix_OV[matrix_OV != 0] <- matrix_OV[matrix_OV != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_OV)
# Reorder the rows of the data frame based on their row sums
matrix_OV <- matrix_OV[order(row_sums, decreasing = FALSE), ]
heatmap_OV <- heatmap_survival(matrix_OV)

#Creating two heatmaps on top of each other
heatmap_OV$fill_color <- ifelse(heatmap_OV$gene.omic == "SCORE" & heatmap_OV$value == 0.2, "lightblue", ifelse(heatmap_OV$gene.omic == "SCORE" & heatmap_OV$value == 0.5, "orange", ifelse(heatmap_OV$gene.omic == "SCORE" & heatmap_OV$value == 0.8, "darkred", scales::col_numeric(c("white", "red"), domain = c(min(heatmap_OV$value), max(heatmap_OV$value)))(heatmap_OV$value))))

p1 <- ggplot(heatmap_OV[heatmap_OV$gene.omic == "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())

p2 <- ggplot(heatmap_OV[heatmap_OV$gene.omic != "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_text(size = 2, angle = 90))

pdf("heatmap_OV.pdf")
plot_grid(p1,p2,nrow=2,align='v',axis='tb',rel_heights=c(1/20,19/20))
dev.off()

### 10.4 UCEC

# Increase all non-zero values by 1 to have better colors
matrix_UCEC[matrix_UCEC != 0] <- matrix_UCEC[matrix_UCEC != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_UCEC)
# Reorder the rows of the data frame based on their row sums
matrix_UCEC <- matrix_UCEC[order(row_sums, decreasing = FALSE), ]
heatmap_UCEC <- heatmap_survival(matrix_UCEC)

#Creating two heatmaps on top of each other
heatmap_UCEC$fill_color <- ifelse(heatmap_UCEC$gene.omic == "SCORE" & heatmap_UCEC$value == 0.2, "lightblue", ifelse(heatmap_UCEC$gene.omic == "SCORE" & heatmap_UCEC$value == 0.5, "orange", ifelse(heatmap_UCEC$gene.omic == "SCORE" & heatmap_UCEC$value == 0.8, "darkred", scales::col_numeric(c("white", "red"), domain = c(min(heatmap_UCEC$value), max(heatmap_UCEC$value)))(heatmap_UCEC$value))))

p1 <- ggplot(heatmap_UCEC[heatmap_UCEC$gene.omic == "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 8), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())

p2 <- ggplot(heatmap_UCEC[heatmap_UCEC$gene.omic != "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 3, angle = 90))

pdf("heatmap_UCEC.pdf")
plot_grid(p1,p2,nrow=2,align='v',axis='tb',rel_heights=c(1/20,19/20))
dev.off()

### 10.5 UCS

# Increase all non-zero values by 1 to have better colors
matrix_UCS[matrix_UCS != 0] <- matrix_UCS[matrix_UCS != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_UCS)
# Reorder the rows of the data frame based on their row sums
matrix_UCS <- matrix_UCS[order(row_sums, decreasing = FALSE), ]
heatmap_UCS <- heatmap_survival(matrix_UCS)

#Creating two heatmaps on top of each other
heatmap_UCS$fill_color <- ifelse(heatmap_UCS$gene.omic == "SCORE" & heatmap_UCS$value == 0.2, "lightblue", ifelse(heatmap_UCS$gene.omic == "SCORE" & heatmap_UCS$value == 0.5, "orange", ifelse(heatmap_UCS$gene.omic == "SCORE" & heatmap_UCS$value == 0.8, "darkred", scales::col_numeric(c("white", "red"), domain = c(min(heatmap_UCS$value), max(heatmap_UCS$value)))(heatmap_UCS$value))))

p1 <- ggplot(heatmap_UCS[heatmap_UCS$gene.omic == "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 12.2), axis.text.x = element_blank(), axis.title.x = element_blank(), axis.title.y = element_blank(), axis.line.x = element_blank(), axis.ticks.x = element_blank())

p2 <- ggplot(heatmap_UCS[heatmap_UCS$gene.omic != "SCORE", ], aes(x = patient, y = gene.omic, fill = fill_color)) +
  geom_tile() +
  scale_fill_identity() +
  theme(axis.text.y = element_text(size = 5), axis.text.x = element_text(size = 5, angle = 90))

pdf("heatmap_UCS.pdf")
plot_grid(p1,p2,nrow=2,align='v',axis='tb',rel_heights=c(1/20,19/20))
dev.off()

## 11. Creating the Kaplan_Meyer plot

KM_fun <- function(df, survival){
  df[df$gene.omic == "SCORE" & df$value == 0.2, "value"] <- "Low"
  df[df$gene.omic == "SCORE" & df$value == 0.5, "value"] <- "Intermediate"
  df[df$gene.omic == "SCORE" & df$value == 0.8, "value"] <- "High"
  df <- df[df$gene.omic == "SCORE", , drop = FALSE]
  rownames(df) <- df$patient
  df <- df[order(row.names(df)), , drop = FALSE]
  survival <- survival[order(row.names(survival)), , drop = FALSE]
  survival$Risk <- df$value
  survival <- survival[order(survival$Risk), , drop = FALSE]
  return(survival)
}

### 11.1 BRCA

survival_BRCA <- KM_fun(heatmap_survival(matrix_BRCA), survival_BRCA)
kmcurve_BRCA <- survfit(Surv(survival_BRCA$days, survival_BRCA$status) ~ Risk, data = survival_BRCA)

pdf("KM_BRCA.pdf")
ggsurvplot(kmcurve_BRCA, xlim = c(0,9000), break.x.by = 1000, pval = TRUE, risk.table = TRUE, palette = c("darkred", "orange", "lightblue"), legend = "right", legend.labs = c("High Risk", "Mid Risk", "Low Risk"))
dev.off()

### 11.2 CESC

survival_CESC <- KM_fun(heatmap_survival(matrix_CESC), survival_CESC)
kmcurve_CESC <- survfit(Surv(survival_CESC$days, survival_CESC$status) ~ Risk, data = survival_CESC)

pdf("KM_CESC.pdf")
ggsurvplot(kmcurve_CESC, xlim = c(0,7000), break.x.by = 1000, pval = TRUE, risk.table = TRUE, palette = c("darkred", "orange", "lightblue"), legend = "right", legend.labs = c("High Risk", "Mid Risk", "Low Risk"))
dev.off()

### 11.3 OV

survival_OV <- KM_fun(heatmap_survival(matrix_OV), survival_OV)
kmcurve_OV <- survfit(Surv(survival_OV$days, survival_OV$status) ~ Risk, data = survival_OV)

pdf("KM_OV.pdf")
ggsurvplot(kmcurve_OV, xlim = c(0,6000), break.x.by = 1000, pval = TRUE, risk.table = TRUE, palette = c("darkred", "orange", "lightblue"), legend = "right", legend.labs = c("High Risk", "Mid Risk", "Low Risk"))
dev.off()

### 11.4 UCEC

survival_UCEC <- KM_fun(heatmap_survival(matrix_UCEC), survival_UCEC)
kmcurve_UCEC <- survfit(Surv(survival_UCEC$days, survival_UCEC$status) ~ Risk, data = survival_UCEC)

pdf("KM_UCEC.pdf")
ggsurvplot(kmcurve_UCEC, xlim = c(0,4000), break.x.by = 1000, pval = TRUE, risk.table = TRUE, palette = c("darkred", "orange", "lightblue"), legend = "right", legend.labs = c("High Risk", "Mid Risk", "Low Risk"))
dev.off()

### 11.5 UCS

survival_UCS <- KM_fun(heatmap_survival(matrix_UCS), survival_UCS)
kmcurve_UCS <- survfit(Surv(survival_UCS$days, survival_UCS$status) ~ Risk, data = survival_UCS)

pdf("KM_UCS.pdf")
ggsurvplot(kmcurve_UCS, xlim = c(0,5000), break.x.by = 1000, pval = TRUE, risk.table = TRUE, palette = c("darkred", "orange", "lightblue"), legend = "right", legend.labs = c("High Risk", "Mid Risk", "Low Risk"))
dev.off()

## 12. Seperating the testing and training

### 12.1 BRCA

#### 12.1.1 Expression

col_BRCA <- sort(sample(1:dim(matrix_BRCA_exp)[2], round((0.1 * dim(matrix_BRCA_exp)[2])), replace = FALSE))
Test_BRCA_EXP <- matrix_BRCA_exp[, col_BRCA]
matrix_BRCA_exp <- matrix_BRCA_exp[, -col_BRCA]

#### 12.1.2 Methylation

Test_BRCA_MET <- matrix_BRCA_met[, col_BRCA]
matrix_BRCA_met <- matrix_BRCA_met[, -col_BRCA]

#### 12.1.3 CNV

Test_BRCA_CNV <- matrix_BRCA_cnv[, col_BRCA]
matrix_BRCA_cnv <- matrix_BRCA_cnv[, -col_BRCA]

#### 12.1.4 Mutataion

Test_BRCA_MUT <- matrix_BRCA_mut[, col_BRCA]
matrix_BRCA_mut <- matrix_BRCA_mut[, -col_BRCA]

#### 12.1.5 Survival

Test_BRCA_survival <- survival_BRCA[rownames(survival_BRCA) %in% colnames(Test_BRCA_EXP), ]
survival_BRCA <- survival_BRCA[!rownames(survival_BRCA) %in% colnames(Test_BRCA_EXP), ]

### 12.2 CESC

#### 12.2.1 Expression

col_CESC <- sort(sample(1:dim(matrix_CESC_exp)[2], round((0.1 * dim(matrix_CESC_exp)[2])), replace = FALSE))
Test_CESC_EXP <- matrix_CESC_exp[, col_CESC]
matrix_CESC_exp <- matrix_CESC_exp[, -col_CESC]

#### 12.2.2 Methylation

Test_CESC_MET <- matrix_CESC_met[, col_CESC]
matrix_CESC_met <- matrix_CESC_met[, -col_CESC]

#### 12.2.3 CNV

Test_CESC_CNV <- matrix_CESC_cnv[, col_CESC]
matrix_CESC_cnv <- matrix_CESC_cnv[, -col_CESC]

#### 12.2.4 Mutataion

Test_CESC_MUT <- matrix_CESC_mut[, col_CESC]
matrix_CESC_mut <- matrix_CESC_mut[, -col_CESC]

#### 12.2.5 Survival

Test_CESC_survival <- survival_CESC[rownames(survival_CESC) %in% colnames(Test_CESC_EXP), ]
survival_CESC <- survival_CESC[!rownames(survival_CESC) %in% colnames(Test_CESC_EXP), ]
Test_CESC_survival

### 12.3 OV

#### 12.3.1 Expression

col_OV <- sort(sample(1:dim(matrix_OV_exp)[2], round((0.1 * dim(matrix_OV_exp)[2])), replace = FALSE))
Test_OV_EXP <- matrix_OV_exp[, col_OV]
matrix_OV_exp <- matrix_OV_exp[, -col_OV]

#### 12.3.2 Methylation

Test_OV_MET <- matrix_OV_met[, col_OV]
matrix_OV_met <- matrix_OV_met[, -col_OV]

#### 12.3.3 CNV

Test_OV_CNV <- matrix_OV_cnv[, col_OV]
matrix_OV_cnv <- matrix_OV_cnv[, -col_OV]

#### 12.3.4 Mutataion

Test_OV_MUT <- matrix_OV_mut[, col_OV]
matrix_OV_mut <- matrix_OV_mut[, -col_OV]

#### 12.3.5 Survival

Test_OV_survival <- survival_OV[rownames(survival_OV) %in% colnames(Test_OV_EXP), ]
survival_OV <- survival_OV[!rownames(survival_OV) %in% colnames(Test_OV_EXP), ]

### 12.4 UCEC

#### 12.4.1 Expression

col_UCEC <- sort(sample(1:dim(matrix_UCEC_exp)[2], round((0.1 * dim(matrix_UCEC_exp)[2])), replace = FALSE))
Test_UCEC_EXP <- matrix_UCEC_exp[, col_UCEC]
matrix_UCEC_exp <- matrix_UCEC_exp[, -col_UCEC]

#### 12.4.2 Methylation

Test_UCEC_MET <- matrix_UCEC_met[, col_UCEC]
matrix_UCEC_met <- matrix_UCEC_met[, -col_UCEC]

#### 12.4.3 CNV

Test_UCEC_CNV <- matrix_UCEC_cnv[, col_UCEC]
matrix_UCEC_cnv <- matrix_UCEC_cnv[, -col_UCEC]

#### 12.4.4 Mutation

Test_UCEC_MUT <- matrix_UCEC_mut[, col_UCEC]
matrix_UCEC_mut <- matrix_UCEC_mut[, -col_UCEC]

#### 12.4.5 Survival

Test_UCEC_survival <- survival_UCEC[rownames(survival_UCEC) %in% colnames(Test_UCEC_EXP), ]
survival_UCEC <- survival_UCEC[!rownames(survival_UCEC) %in% colnames(Test_UCEC_EXP), ]

### 12.5 UCS

#### 12.5.1 Expression

col_UCS <- sort(sample(1:dim(matrix_UCS_exp)[2], round((0.1 * dim(matrix_UCS_exp)[2])), replace = FALSE))
Test_UCS_EXP <- matrix_UCS_exp[, col_UCS]
matrix_UCS_exp <- matrix_UCS_exp[, -col_UCS]

#### 12.5.2 Methylation

Test_UCS_MET <- matrix_UCS_met[, col_UCS]
matrix_UCS_met <- matrix_UCS_met[, -col_UCS]

#### 12.5.3 CNV

Test_UCS_CNV <- matrix_UCS_cnv[, col_UCS]
matrix_UCS_cnv <- matrix_UCS_cnv[, -col_UCS]

#### 12.5.4 Mutation

Test_UCS_MUT <- matrix_UCS_mut[, col_UCS]
matrix_UCS_mut <- matrix_UCS_mut[, -col_UCS]

#### 12.5.5 Survival

Test_UCS_survival <- survival_UCS[rownames(survival_UCS) %in% colnames(Test_UCS_EXP), ]
survival_UCS <- survival_UCS[!rownames(survival_UCS) %in% colnames(Test_UCS_EXP), ]

## 13. Training the data 

#A function for expression and methylation
exp_met_matrix_train_test <- function(survival, matrix, Test, l){
  penalized <- penalized::penalized(Surv(survival$days, survival$status), penalized = as.data.frame(t(matrix)), unpenalized = ~0, lambda1 = l, model = "cox")
  HR <- as.data.frame(coefficients(penalized))
  # Here after training the data we added the testing data to the trained data
  matrix <- merge(matrix, Test, by = 0, all = TRUE)
  rownames(matrix) <- matrix$Row.names
  matrix$Row.names <- NULL
  matrix <- matrix[rownames(matrix) %in% rownames(HR), ]
  for (gene in rownames(matrix)){
    median <- median(as.numeric(matrix[gene,]))
    values <- as.numeric(matrix[gene,])
    max <- (max(values) + 1)
    if (HR[gene,] < 0){
      values[values >= median] <- max
      values[values < median] <- abs(as.numeric(HR[gene,]))
      values[values == max] <- 0
    }  
    else if(HR[gene,] > 0){
      values[values > median] <- max
      values[values <= median] <-  0
      values[values == max] <- abs(as.numeric(HR[gene,]))
    }
    matrix[gene,] <- values
  }
  return(matrix)
}

#A function for CNV
cnv_matrix_train_test <- function(survival, matrix, Test, l){
  penalized <- penalized::penalized(Surv(survival$days, survival$status), penalized = as.data.frame(t(matrix)), unpenalized = ~0, lambda1 = l, model = "cox")
  HR <- as.data.frame(coefficients(penalized))
  # Here after training the data we added the testing data to the trained data
  matrix <- merge(matrix, Test, by = 0, all = TRUE)
  rownames(matrix) <- matrix$Row.names
  matrix$Row.names <- NULL
  matrix <- matrix[rownames(matrix) %in% rownames(HR), ]
  for (gene in rownames(matrix)){
    values <- as.numeric(matrix[gene,])
    if (HR[gene,] < 0){
      values[values >= 0] <- 0
      values[values < 0] <- abs(as.numeric(HR[gene,]))
    }
    else if(HR[gene,] > 0){
      values[values > 0] <- abs(as.numeric(HR[gene,]))
      values[values <= 0] <- 0
    }
    matrix[gene,] <- values
  }
  return(matrix)
}

#A function for mutation
mutation_matrix_train_test <- function(survival, matrix, Test, l){
  penalized <- penalized::penalized(Surv(survival$days, survival$status), penalized = as.data.frame(t(matrix)), unpenalized = ~0, lambda1 = l, model = "cox")
  HR <- as.data.frame(coefficients(penalized))
  # Here after training the data we added the testing data to the trained data
  matrix <- merge(matrix, Test, by = 0, all = TRUE)
  rownames(matrix) <- matrix$Row.names
  matrix$Row.names <- NULL
  matrix <- matrix[rownames(matrix) %in% rownames(HR), ]
  for (gene in rownames(matrix)){
    values <- as.numeric(matrix[gene,])
    if (HR[gene,] < 0){
      values[values == 1] <- 0
    }
    else if(HR[gene,] > 0){
      values[values == 1] <- abs(as.numeric(HR[gene,]))
    }
    matrix[gene,] <- values
  }
  return(matrix)
}

### 13.1 BRCA

#### 13.1.1 Expression

BRCA_exp_final <- exp_met_matrix_train_test(survival_BRCA, matrix_BRCA_exp, Test_BRCA_EXP, 3.4)

#### 13.1.2 Methylation

BRCA_met_final <- exp_met_matrix_train_test(survival_BRCA, matrix_BRCA_met, Test_BRCA_MET, 3.4)

#### 13.1.3 CNV

BRCA_cnv_final <- cnv_matrix_train_test(survival_BRCA, matrix_BRCA_cnv, Test_BRCA_CNV, 3.4)

#### 13.1.4 Mutation

BRCA_mut_final <- mutation_matrix_train_test(survival_BRCA, matrix_BRCA_mut, Test_BRCA_MUT, 3.4)

matrix_BRCA_Training <- do.call("rbind", list(BRCA_exp_final, BRCA_cnv_final, BRCA_met_final, BRCA_mut_final))
save(matrix_BRCA_Training, file = "./matrix_BRCA_Training.RData")

matrix_BRCA_Testing <- do.call("rbind", list(Test_BRCA_EXP, Test_BRCA_CNV, Test_BRCA_MET, Test_BRCA_MUT))
save(matrix_BRCA_Testing, file = "./matrix_BRCA_Testing.RData")

### 13.2 CESC

#### 13.2.1 Expression

CESC_exp_final <- exp_met_matrix_train_test(survival_CESC, matrix_CESC_exp, Test_CESC_EXP, 1.4)

#### 13.2.2 Methylation

CESC_met_final <- exp_met_matrix_train_test(survival_CESC, matrix_CESC_met, Test_CESC_MET, 1.4)

#### 13.2.3 CNV

CESC_cnv_final <- cnv_matrix_train_test(survival_CESC, matrix_CESC_cnv, Test_CESC_CNV, 1.4)

#### 13.2.4 Mutataion

CESC_mut_final <- mutation_matrix_train_test(survival_CESC, matrix_CESC_mut, Test_CESC_MUT, 1.4)

matrix_CESC_Training <- do.call("rbind", list(CESC_exp_final, CESC_cnv_final, CESC_met_final, CESC_mut_final))
save(matrix_CESC_Training, file = "./matrix_CESC_Training.RData")

matrix_CESC_Testing <- do.call("rbind", list(Test_CESC_EXP, Test_CESC_CNV, Test_CESC_MET, Test_CESC_MUT))
save(matrix_CESC_Testing, file = "./matrix_CESC_Testing.RData")

### 13.3 OV

#### 13.3.1 Expression

OV_exp_final <- exp_met_matrix_train_test(survival_OV, matrix_OV_exp, Test_OV_EXP,3.3)

#### 13.3.2 Methylation

OV_met_final <- exp_met_matrix_train_test(survival_OV, matrix_OV_met, Test_OV_MET, 3.3)

#### 13.3.3 CNV

OV_cnv_final <- cnv_matrix_train_test(survival_OV, matrix_OV_cnv, Test_OV_CNV, 3.3)

#### 13.3.4 Mutataion

OV_mut_final <- mutation_matrix_train_test(survival_OV, matrix_OV_mut, Test_OV_MUT, 3.3)

matrix_OV_Training <- do.call("rbind", list(OV_exp_final, OV_cnv_final, OV_met_final, OV_mut_final))
save(matrix_OV_Training, file = "./matrix_OV_Training.RData")

matrix_OV_Testing <- do.call("rbind", list(Test_OV_EXP, Test_OV_CNV, Test_OV_MET, Test_OV_MUT))
save(matrix_OV_Testing, file = "./matrix_OV_Testing.RData")

### 13.4 UCEC

#### 13.4.1 Expression

UCEC_exp_final <- exp_met_matrix_train_test(survival_UCEC, matrix_UCEC_exp, Test_UCEC_EXP, 3)

#### 13.4.2 Methylation

UCEC_met_final <- exp_met_matrix_train_test(survival_UCEC, matrix_UCEC_met, Test_UCEC_MET, 3)

#### 13.4.3 CNV

UCEC_cnv_final <- cnv_matrix_train_test(survival_UCEC, matrix_UCEC_cnv, Test_UCEC_CNV, 3)

#### 13.4.4 Mutataion

UCEC_mut_final <- mutation_matrix_train_test(survival_UCEC, matrix_UCEC_mut, Test_UCEC_MUT, 3)

matrix_UCEC_Training <- do.call("rbind", list(UCEC_exp_final, UCEC_cnv_final, UCEC_met_final, UCEC_mut_final))
save(matrix_UCEC_Training, file = "./matrix_UCEC_Training.RData")

matrix_UCEC_Testing <- do.call("rbind", list(Test_UCEC_EXP, Test_UCEC_CNV, Test_UCEC_MET, Test_UCEC_MUT))
save(matrix_UCEC_Testing, file = "./matrix_UCEC_Testing.RData")

### 13.5 UCS

#### 13.5.1 Expression

UCS_exp_final <- exp_met_matrix_train_test(survival_UCS, matrix_UCS_exp, Test_UCS_EXP, 1.5)

#### 13.5.2 Methylation

UCS_met_final <- exp_met_matrix_train_test(survival_UCS, matrix_UCS_met, Test_UCS_MET, 1.5)

#### 13.5.3 CNV

UCS_cnv_final <- cnv_matrix_train_test(survival_UCS, matrix_UCS_cnv, Test_UCS_CNV, 1.5)

#### 13.5.4 Mutataion

UCS_mut_final <- mutation_matrix_train_test(survival_UCS, matrix_UCS_mut, Test_UCS_MUT, 1.5)

matrix_UCS_Training <- do.call("rbind", list(UCS_exp_final, UCS_cnv_final, UCS_met_final, UCS_mut_final))
save(matrix_UCS_Training, file = "./matrix_UCS_Training.RData")

matrix_UCS_Testing <- do.call("rbind", list(Test_UCS_EXP, Test_UCS_CNV, Test_UCS_MET, Test_UCS_MUT))
save(matrix_UCS_Testing, file = "./matrix_UCS_Testing.RData")

## 14. Testing the data

heatmap_survival_train_test <- function(matrix){
  #changing the gene id to symbol
  genes <- list()
  for (gene in rownames(matrix)){
    b <- strsplit(gene, ":")
    c <- b[[1]][2]
    d <- strsplit(c, "[.]")
    genes <- c(genes, d[[1]][1])
  }
  genes <- as.character(genes)
  symbol <- AnnotationDbi::select(org.Hs.eg.db, keys= genes, keytype = "ENTREZID" ,columns = "SYMBOL")
  matrix$symbol <- symbol$SYMBOL
  matrix <- matrix %>% relocate(symbol)
  n <- 0
  for (gene in rownames(matrix)){
    n <- n+1
    d <- strsplit(gene, "[.]")
    matrix$symbol[n] <- paste0(matrix$symbol[n], ".", d[[1]][2])
  }
  rownames(matrix) <- matrix$symbol
  matrix$symbol <- NULL
  #making the total survival score
  patients_matrix <- colnames(matrix)
  for (patient in patients_matrix){
    matrix["SCORE", patient] <- 0
    matrix["SCORE", patient] <- sum(matrix[, patient])
  }
  matrix <- as.data.frame(t(matrix))
  matrix <- matrix[order(matrix$SCORE),]
  #rescaling the scores value
  q <- as.list(quantile(matrix$SCORE))
  matrix[matrix$SCORE <= q[[2]], "SCORE"] <- 0.2
  matrix[matrix$SCORE > q[[2]] & matrix$SCORE <= q[[4]], "SCORE"] <- 0.5
  matrix[matrix$SCORE > q[[4]] & matrix$SCORE <= q[[5]], "SCORE"] <- 0.8
  return(matrix)
  }

### 14.1 BRCA

#### 14.1.1 Tested data

# Increase all non-zero values by 1 to have better colors
matrix_BRCA_Training[matrix_BRCA_Training != 0] <- matrix_BRCA_Training[matrix_BRCA_Training != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_BRCA_Training)
# Reorder the rows of the data frame based on their row sums
matrix_BRCA_Training <- matrix_BRCA_Training[order(row_sums, decreasing = FALSE), ]
heatmap_BRCA_train_test <- heatmap_survival_train_test(matrix_BRCA_Training)
heatmap_BRCA_train_test[,"SCORE", drop=F]

BRCA_Testing <- heatmap_BRCA_train_test[rownames(heatmap_BRCA_train_test) %in% colnames(matrix_BRCA_Testing), "SCORE", drop = F]
# Change values of cells
BRCA_Testing$SCORE[BRCA_Testing$SCORE == 0.2] <- "Low"
BRCA_Testing$SCORE[BRCA_Testing$SCORE == 0.5] <- "Mid"
BRCA_Testing$SCORE[BRCA_Testing$SCORE == 0.8] <- "High"
save(BRCA_Testing, file = "./BRCA_Testing.RData")

#### 14.1.2 Old data

# Increase all non-zero values by 1 to have better colors
matrix_BRCA[matrix_BRCA != 0] <- matrix_BRCA[matrix_BRCA != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_BRCA)
# Reorder the rows of the data frame based on their row sums
matrix_BRCA <- matrix_BRCA[order(row_sums, decreasing = FALSE), ]
heatmap_BRCA <- heatmap_survival_train_test(matrix_BRCA)
heatmap_BRCA[,"SCORE", drop=F]

heatmap_BRCA <- heatmap_BRCA[rownames(heatmap_BRCA) %in% colnames(matrix_BRCA_Testing), "SCORE", drop = F]
# Change values of cells
heatmap_BRCA$SCORE[heatmap_BRCA$SCORE == 0.2] <- "Low"
heatmap_BRCA$SCORE[heatmap_BRCA$SCORE == 0.5] <- "Mid"
heatmap_BRCA$SCORE[heatmap_BRCA$SCORE == 0.8] <- "High"

# Calculate percentage of similarity
similarity_BRCA <- sum(BRCA_Testing$SCORE == heatmap_BRCA[match(rownames(BRCA_Testing), rownames(heatmap_BRCA)), "SCORE"]) / nrow(BRCA_Testing) * 100

### 14.2 CESC

#### 14.2.1 Tested data

# Increase all non-zero values by 1 to have better colors
matrix_CESC_Training[matrix_CESC_Training != 0] <- matrix_CESC_Training[matrix_CESC_Training != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_CESC_Training)
# Reorder the rows of the data frame based on their row sums
matrix_CESC_Training <- matrix_CESC_Training[order(row_sums, decreasing = FALSE), ]
heatmap_CESC_train_test <- heatmap_survival_train_test(matrix_CESC_Training)
heatmap_CESC_train_test[,"SCORE", drop=F]

CESC_Testing <- heatmap_CESC_train_test[rownames(heatmap_CESC_train_test) %in% colnames(matrix_CESC_Testing), "SCORE", drop = F]
# Change values of cells
CESC_Testing$SCORE[CESC_Testing$SCORE == 0.2] <- "Low"
CESC_Testing$SCORE[CESC_Testing$SCORE == 0.5] <- "Mid"
CESC_Testing$SCORE[CESC_Testing$SCORE == 0.8] <- "High"
save(CESC_Testing, file = "./CESC_Testing.RData")

#### 14.2.2 Old data

# Increase all non-zero values by 1 to have better colors
matrix_CESC[matrix_CESC != 0] <- matrix_CESC[matrix_CESC != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_CESC)
# Reorder the rows of the data frame based on their row sums
matrix_CESC <- matrix_CESC[order(row_sums, decreasing = FALSE), ]
heatmap_CESC <- heatmap_survival_train_test(matrix_CESC)
heatmap_CESC[,"SCORE", drop=F]

heatmap_CESC <- heatmap_CESC[rownames(heatmap_CESC) %in% colnames(matrix_CESC_Testing), "SCORE", drop = F]
# Change values of cells
heatmap_CESC$SCORE[heatmap_CESC$SCORE == 0.2] <- "Low"
heatmap_CESC$SCORE[heatmap_CESC$SCORE == 0.5] <- "Mid"
heatmap_CESC$SCORE[heatmap_CESC$SCORE == 0.8] <- "High"

# Calculate percentage of similarity
similarity_CESC <- sum(CESC_Testing$SCORE == heatmap_CESC[match(rownames(CESC_Testing), rownames(heatmap_CESC)), "SCORE"]) / nrow(CESC_Testing) * 100

### 14.3 OV

#### 14.3.1 Tested data

# Increase all non-zero values by 1 to have better colors
matrix_OV_Training[matrix_OV_Training != 0] <- matrix_OV_Training[matrix_OV_Training != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_OV_Training)
# Reorder the rows of the data frame based on their row sums
matrix_OV_Training <- matrix_OV_Training[order(row_sums, decreasing = FALSE), ]
heatmap_OV_train_test <- heatmap_survival_train_test(matrix_OV_Training)
heatmap_OV_train_test[,"SCORE", drop=F]

OV_Testing <- heatmap_OV_train_test[rownames(heatmap_OV_train_test) %in% colnames(matrix_OV_Testing), "SCORE", drop = F]
# Change values of cells
OV_Testing$SCORE[OV_Testing$SCORE == 0.2] <- "Low"
OV_Testing$SCORE[OV_Testing$SCORE == 0.5] <- "Mid"
OV_Testing$SCORE[OV_Testing$SCORE == 0.8] <- "High"
save(OV_Testing, file = "./OV_Testing.RData")

#### 14.1.2 Old data

# Increase all non-zero values by 1 to have better colors
matrix_OV[matrix_OV != 0] <- matrix_OV[matrix_OV != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_OV)
# Reorder the rows of the data frame based on their row sums
matrix_OV <- matrix_OV[order(row_sums, decreasing = FALSE), ]
heatmap_OV <- heatmap_survival_train_test(matrix_OV)
heatmap_OV[,"SCORE", drop=F]

heatmap_OV <- heatmap_OV[rownames(heatmap_OV) %in% colnames(matrix_OV_Testing), "SCORE", drop = F]
# Change values of cells
heatmap_OV$SCORE[heatmap_OV$SCORE == 0.2] <- "Low"
heatmap_OV$SCORE[heatmap_OV$SCORE == 0.5] <- "Mid"
heatmap_OV$SCORE[heatmap_OV$SCORE == 0.8] <- "High"

# Calculate percentage of similarity
similarity_OV <- sum(OV_Testing$SCORE == heatmap_OV[match(rownames(OV_Testing), rownames(heatmap_OV)), "SCORE"]) / nrow(OV_Testing) * 100

### 14.4 UCEC

#### 14.4.1 Tested data

# Increase all non-zero values by 1 to have better colors
matrix_UCEC_Training[matrix_UCEC_Training != 0] <- matrix_UCEC_Training[matrix_UCEC_Training != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_UCEC_Training)
# Reorder the rows of the data frame based on their row sums
matrix_UCEC_Training <- matrix_UCEC_Training[order(row_sums, decreasing = FALSE), ]
heatmap_UCEC_train_test <- heatmap_survival_train_test(matrix_UCEC_Training)
heatmap_UCEC_train_test[,"SCORE", drop=F]

UCEC_Testing <- heatmap_UCEC_train_test[rownames(heatmap_UCEC_train_test) %in% colnames(matrix_UCEC_Testing), "SCORE", drop = F]
# Change values of cells
UCEC_Testing$SCORE[UCEC_Testing$SCORE == 0.2] <- "Low"
UCEC_Testing$SCORE[UCEC_Testing$SCORE == 0.5] <- "Mid"
UCEC_Testing$SCORE[UCEC_Testing$SCORE == 0.8] <- "High"
save(UCEC_Testing, file = "./UCEC_Testing.RData")

#### 14.4.2 Old data

# Increase all non-zero values by 1 to have better colors
matrix_UCEC[matrix_UCEC != 0] <- matrix_UCEC[matrix_UCEC != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_UCEC)
# Reorder the rows of the data frame based on their row sums
matrix_UCEC <- matrix_UCEC[order(row_sums, decreasing = FALSE), ]
heatmap_UCEC <- heatmap_survival_train_test(matrix_UCEC)
heatmap_UCEC[,"SCORE", drop=F]

heatmap_UCEC <- heatmap_UCEC[rownames(heatmap_UCEC) %in% colnames(matrix_UCEC_Testing), "SCORE", drop = F]
# Change values of cells
heatmap_UCEC$SCORE[heatmap_UCEC$SCORE == 0.2] <- "Low"
heatmap_UCEC$SCORE[heatmap_UCEC$SCORE == 0.5] <- "Mid"
heatmap_UCEC$SCORE[heatmap_UCEC$SCORE == 0.8] <- "High"

# Calculate percentage of similarity
similarity_UCEC <- sum(UCEC_Testing$SCORE == heatmap_UCEC[match(rownames(UCEC_Testing), rownames(heatmap_UCEC)), "SCORE"]) / nrow(UCEC_Testing) * 100

### 14.5 UCS

#### 14.5.1 Tested data

# Increase all non-zero values by 1 to have better colors
matrix_UCS_Training[matrix_UCS_Training != 0] <- matrix_UCS_Training[matrix_UCS_Training != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_UCS_Training)
# Reorder the rows of the data frame based on their row sums
matrix_UCS_Training <- matrix_UCS_Training[order(row_sums, decreasing = FALSE), ]
heatmap_UCS_train_test <- heatmap_survival_train_test(matrix_UCS_Training)
heatmap_UCS_train_test[,"SCORE", drop=F]

UCS_Testing <- heatmap_UCS_train_test[rownames(heatmap_UCS_train_test) %in% colnames(matrix_UCS_Testing), "SCORE", drop = F]
# Change values of cells
UCS_Testing$SCORE[UCS_Testing$SCORE == 0.2] <- "Low"
UCS_Testing$SCORE[UCS_Testing$SCORE == 0.5] <- "Mid"
UCS_Testing$SCORE[UCS_Testing$SCORE == 0.8] <- "High"
save(UCS_Testing, file = "./UCS_Testing.RData")

#### 14.5.2 Old data

# Increase all non-zero values by 1 to have better colors
matrix_UCS[matrix_UCS != 0] <- matrix_UCS[matrix_UCS != 0] + 0.1
# Calculate the row sums
row_sums <- rowSums(matrix_UCS)
# Reorder the rows of the data frame based on their row sums
matrix_UCS <- matrix_UCS[order(row_sums, decreasing = FALSE), ]
heatmap_UCS <- heatmap_survival_train_test(matrix_UCS)
heatmap_UCS[,"SCORE", drop=F]

heatmap_UCS <- heatmap_UCS[rownames(heatmap_UCS) %in% colnames(matrix_UCS_Testing), "SCORE", drop = F]
# Change values of cells
heatmap_UCS$SCORE[heatmap_UCS$SCORE == 0.2] <- "Low"
heatmap_UCS$SCORE[heatmap_UCS$SCORE == 0.5] <- "Mid"
heatmap_UCS$SCORE[heatmap_UCS$SCORE == 0.8] <- "High"

# Calculate percentage of similarity
similarity_UCS <- sum(UCS_Testing$SCORE == heatmap_UCS[match(rownames(UCS_Testing), rownames(heatmap_UCS)), "SCORE"]) / nrow(UCS_Testing) * 100

