# Load Environment Libraries
library(edgeR)
library(limma)
library(biomaRt)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(tidyverse)
library(plotly)
library(RColorBrewer)
library(gameofthrones)
library(gprofiler2)
library(readxl)
library(matrixStats)
library(cowplot)
library(tibble)
library(org.Hs.eg.db)

# Set Working Directory
setwd("C:\\Users\\DELL\\OneDrive\\Documents\\RNA_SEQ_CRC")

# Load Raw Counts Matrix
rna_data <- read_csv("counts.csv")
# Clean and Standardize Sample Columns
rna_data$Normal1 <- rna_data$S1_Normal_sorted.bam
rna_data$Normal2 <- rna_data$S2_Normal_sorted.bam
rna_data$Normal3 <- rna_data$S3_Normal_sorted.bam
rna_data$Normal4 <- rna_data$S4_Normal_sorted.bam
rna_data$Normal5 <- rna_data$S5_Normal_sorted.bam
rna_data$Tumor1  <- rna_data$S1_Tumor_sorted.bam
rna_data$Tumor2  <- rna_data$S2_Tumor_sorted.bam
rna_data$Tumor3  <- rna_data$S3_Tumor_sorted.bam
rna_data$Tumor4  <- rna_data$S4_Tumor_sorted.bam
rna_data$Tumor5  <- rna_data$S5_Tumor_sorted.bam

rna_data <- select(rna_data, Geneid, Normal1, Normal2, Normal3, Normal4, Normal5, Tumor1, Tumor2, Tumor3, Tumor4, Tumor5)

# Preparing Count Matrix Dataframe
mycounts <- as.data.frame(rna_data)
names(mycounts)[names(mycounts) == "Geneid"] <- "Gene_ID"
rownames(mycounts) <- mycounts$Gene_ID
mycounts <- mycounts[ , -1]

# PAIRED DESIGN CONFIGURATION: Set up factors tracking both Condition and Patient
group <- factor(c("Normal", "Normal", "Normal", "Normal", "Normal", "Tumor", "Tumor", "Tumor", "Tumor", "Tumor"))
patient <- factor(c("P1", "P2", "P3", "P4", "P5", "P1", "P2", "P3", "P4", "P5"))

targets <- data.frame(row.names = colnames(mycounts), group, patient)
sample_labels <- rownames(targets)

# Instantiate edgeR DGEList Object
myDGElist <- DGEList(mycounts)

#Data processing for GSEA
myDGElist.unfiltered <- myDGElist
log2.cpm.unfiltered <- cpm(myDGElist.unfiltered, log=TRUE)
log2.cpm.unfiltered.df <- as_tibble(log2.cpm.unfiltered, rownames= "Gene_ID")
colnames(log2.cpm.unfiltered.df) <- c("Gene_ID", sample_labels)
log2.cpm.unfiltered.df.pivot <- pivot_longer(log2.cpm.unfiltered.df, cols= -1, names_to = "Samples", values_to = "Expression")

#Normalization
myDGElist.unfiltered.norm <- calcNormFactors(myDGElist.unfiltered, method = "TMM")
log2.cpm.unfiltered.norm <- cpm(myDGElist.unfiltered.norm, log=TRUE)
log2.cpm.unfiltered.norm.df <- as_tibble(log2.cpm.unfiltered.norm, rownames= "Gene_ID")
colnames(log2.cpm.unfiltered.norm.df) <- c("Gene_ID", sample_labels)
log2.cpm.unfiltered.norm.df.pivot <- pivot_longer(log2.cpm.unfiltered.norm.df, cols= -1, names_to = "Samples", values_to = "Expression")

# Paired Model Design Matrix Creation
design <- model.matrix(~0 + group + patient)
colnames(design) <- gsub("group", "", colnames(design)) # Clean names for contrast evaluation

# Voom Transformation to account for Mean-Variance relationships
voom_result_unf <- voom(myDGElist.unfiltered.norm, design)

# Fitting Linear Model
model_fit_unf <- lmFit(voom_result_unf, design)

# Evaluate Contrast Matrix (Tumor vs Normal Baseline)
contrast_matrix_unf <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
linear_model_unf <- contrasts.fit(model_fit_unf, contrast_matrix_unf)
bay_stat_unf <- eBayes(linear_model_unf)

# Extract Complete Stats Table (Referencing explicit contrast name)
my_genes_unf <- topTable(bay_stat_unf, adjust = "BH", coef = "Tumor_vs_Normal", number = Inf, sort.by = "P")

# Annotation Phase: Query Ensembl biomaRt for Human Symbols
Ensemble_id_unf <- gsub("\\..*", "", row.names(my_genes_unf))
my_genes_unf$modified_ensembl <- Ensemble_id_unf

ensemble_unf <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
result_unf <- getBM(mart = ensemble_unf,
                    attributes = c("ensembl_gene_id", "hgnc_symbol"),
                    filters = "ensembl_gene_id",
                    values = Ensemble_id_unf)

# Merge mappings into Dataframe without dropping target columns
ann_degunf <- merge(my_genes_unf, result_unf, by.x = "modified_ensembl", by.y = "ensembl_gene_id", all.x = TRUE)
# Handle any blank/empty strings by converting them to NA
ann_degunf$hgnc_symbol[ann_degunf$hgnc_symbol == ""] <- NA
# replace rows where hgnc_symbol is NA with modified_ensembl values
ann_degunf$hgnc_symbol[is.na(ann_degunf$hgnc_symbol)] <- ann_degunf$modified_ensembl[is.na(ann_degunf$hgnc_symbol)]
#Check for duplicates
dup_geneid <- ann_degunf$hgnc_symbol[duplicated(ann_degunf$hgnc_symbol) | duplicated(ann_degunf$hgnc_symbol, fromLast= TRUE)]
#Removing duplicates
ann_degunf <- ann_degunf[!ann_degunf$hgnc_symbol %in% dup_geneid, ]
# set ann_deg row names to hgnc_symbol
row.names(ann_degunf) <- ann_degunf$hgnc_symbol
# remove modified_ensembl and hgnc_symbol columns
ann_degunf <- ann_degunf[, !colnames(ann_degunf)%in%c("modified_ensembl", "hgnc_symbol")]

#GSEA
data <- ann_degunf %>% dplyr::arrange(desc(logFC))
gene_list <- data$logFC
names(gene_list) <- row.names(data)

gse <- gseGO(geneList = gene_list, ont = "All", keyType = "SYMBOL", nPerm = 1000, minGSSize = 3, maxGSSize = 100, pvalueCutoff = 0.05, verbose = TRUE,
             OrgDb = org.Hs.eg.db, pAdjustMethod = "none")

#To filter by NES
gse@result <- gse@result %>% arrange(desc(NES))
view(summary(gse))
gse_result <- gse@result
write.csv(gse_result, file = "GSEA_result.csv")
#dotplot
png("gse_dotplot.png", width = 500, height = 520)
dotplot(gse, showCategory = 10, split = ".sign") + facet_grid(.~.sign) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 7))
dev.off()

#Plot specific GSE result
png("gse1.png", width = 500, height = 400)
#Find the index, NES, and p.adjust value for the desired gene set
gene_set_index <- which(gse@result$Description == "hypothalamus cell differentiation")
gene_set_NES <- gse@result$NES[gene_set_index]
gene_set_p.adj <- gse@result$p.adjust[gene_set_index]

#Create the plot with NES and p.adjust included in the title
gseaplot2(gse, geneSetID = gene_set_index, title = paste("hypothalamus cell differentiation", "\nNES =", round(gene_set_NES, 2),
                                                         ", p.adjust =", formatC(gene_set_p.adj, format = "e", digits =2)))
dev.off()
png("gse2.png", width = 500, height = 400)
gene_set_index <- which(gse@result$Description == "hypersensitivity")
gene_set_NES <- gse@result$NES[gene_set_index]
gene_set_p.adj <- gse@result$p.adjust[gene_set_index]

#Create the plot with NES and p.adjust included in the title
gseaplot2(gse, geneSetID = gene_set_index, title = paste("hypersensitivity", "\nNES =", round(gene_set_NES, 2),
                                                         ", p.adjust =", formatC(gene_set_p.adj, format = "e", digits =2)))
dev.off()