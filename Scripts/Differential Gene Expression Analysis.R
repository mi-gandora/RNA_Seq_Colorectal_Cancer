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
head(rna_data, 5)

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

# Low-Count Filtering and TMM Normalization
cpm <- cpm(myDGElist)
filt_threshold <- rowSums(cpm > 1) >= 2
myDGElist.filtered <- myDGElist[filt_threshold, ]
myDGElist.filtered.norm <- calcNormFactors(myDGElist.filtered, method = "TMM")

# Paired Model Design Matrix Creation
design <- model.matrix(~0 + group + patient)
colnames(design) <- gsub("group", "", colnames(design)) # Clean names for contrast evaluation

# Voom Transformation to account for Mean-Variance relationships
voom_result <- voom(myDGElist.filtered.norm, design)

# Fitting Linear Model
model_fit <- lmFit(voom_result, design)

# Evaluate Contrast Matrix (Tumor vs Normal Baseline)
contrast_matrix <- makeContrasts(Tumor_vs_Normal = Tumor - Normal, levels = design)
linear_model <- contrasts.fit(model_fit, contrast_matrix)
bay_stat <- eBayes(linear_model)

# Extract Complete Stats Table (Referencing explicit contrast name)
my_genes <- topTable(bay_stat, adjust = "BH", coef = "Tumor_vs_Normal", number = Inf, sort.by = "P")

# Extract Significantly Quantified DEGs
diff_exp_genes <- my_genes[(my_genes$adj.P.Val < 0.05 & my_genes$logFC > 1) | (my_genes$adj.P.Val < 0.05 & my_genes$logFC < -1), ]

# Annotation Phase: Query Ensembl biomaRt for Human Symbols
Ensemble_id <- gsub("\\..*", "", row.names(diff_exp_genes))
diff_exp_genes$modified_ensembl <- Ensemble_id

ensemble <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
result <- getBM(mart = ensemble,
                attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = "ensembl_gene_id",
                values = Ensemble_id)

# Merge mappings into Dataframe without dropping target columns
ann_deg <- merge(diff_exp_genes, result, by.x = "modified_ensembl", by.y = "ensembl_gene_id", all.x = TRUE)
# Handle any blank/empty strings by converting them to NA
ann_deg$hgnc_symbol[ann_deg$hgnc_symbol == ""] <- NA
# replace rows where hgnc_symbol is NA with modified_ensembl values
ann_deg$hgnc_symbol[is.na(ann_deg$hgnc_symbol)] <- ann_deg$modified_ensembl[is.na(ann_deg$hgnc_symbol)]
# Force duplicate gene symbols to be unique
unique_symbols <- make.unique(ann_deg$hgnc_symbol)
# set ann_deg row names to hgnc_symbol
row.names(ann_deg) <- ann_deg$hgnc_symbol
# remove modified_ensembl and hgnc_symbol columns
ann_deg <- ann_deg[, !colnames(ann_deg)%in%c("modified_ensembl", "hgnc_symbol")]
# Make the gene symbol the first column
ann_deg$Gene_Symbol <- rownames(ann_deg)
ann_deg <- ann_deg %>% select(Gene_Symbol, everything())

# Separate Up and Downregulated Genes
upregulated_genes <- ann_deg[ann_deg$adj.P.Val < 0.05 & ann_deg$logFC > 1, ]
downregulated_genes <- ann_deg[ann_deg$adj.P.Val < 0.05 & ann_deg$logFC < -1, ]

# Write Pipeline Exports to Disk
write.csv(ann_deg, file = "total_degs.csv", row.names = FALSE)
write.csv(upregulated_genes, file = "upregulated_genes.csv", row.names = FALSE)
write.csv(downregulated_genes, file = "downregulated_genes.csv", row.names = FALSE)

# Generate Visualization (Volcano Plot)
degs <- my_genes %>% as_tibble(rownames = "Gene_ID")

png("volcanoplots.png", width = 500, height = 500)
ggplot(degs) + aes(y = -log10(adj.P.Val), x = logFC) +
  geom_point(aes(color = ifelse(adj.P.Val < 0.05 & logFC > 1, "Upregulated", 
                                ifelse(adj.P.Val < 0.05 & logFC < -1, "Downregulated", "Insignificant"))), size = 2) +
  geom_hline(yintercept = -log10(0.05), linetype = "longdash", color = "grey", linewidth = 1) +
  geom_vline(xintercept = 1, linetype = "longdash", color = "grey", linewidth = 1) +
  geom_vline(xintercept = -1, linetype = "longdash", color = "grey", linewidth = 1) + 
  labs(title = "Tumor vs Matched Adjacent Normal", subtitle = "Differential Gene Expression Analysis") +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "green", "Insignificant" = "grey")) +
  theme_bw() + 
  theme(legend.title = element_blank())
dev.off()

#Gene Ontology
human_mart <- useMart(biomart = "ensembl", host = "ensembl.org", dataset = "hsapiens_gene_ensembl")
ann_diff <- getBM(values = row.names(ann_deg),mart = human_mart, attributes = c("hgnc_symbol", "entrezgene_id", "description"), filters = "hgnc_symbol")
ann_diff$entrezgene_id <- as.character(ann_diff$entrezgene_id)

#Biological Process Analysis
ora_analysis_bp <- enrichGO(gene = ann_diff$entrezgene_id, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP",
                            pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = FALSE, pool = FALSE)
ora_analysis_bp_df <- as.data.frame(ora_analysis_bp)
ora_analysis_bp_final <- clusterProfiler::simplify(ora_analysis_bp)
write_delim(x = as.data.frame(ora_analysis_bp@result), file = "Biological Process.csv", delim = ",")

#dotplot
png("Bp dotplot.png", width = 460, height = 550)
dotplot(ora_analysis_bp_final, showCategory = 10)
dev.off()

#Cellular Component
ora_analysis_cc <- enrichGO(gene = ann_diff$entrezgene_id, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC",
                            pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = FALSE, pool = FALSE)
ora_analysis_cc_df <- as.data.frame(ora_analysis_cc)
ora_analysis_cc_final <- clusterProfiler::simplify(ora_analysis_cc)
write_delim(x = as.data.frame(ora_analysis_cc@result), file = "Cellular Component.csv", delim = ",")
#dotplot
png("CC dotplot.png", width = 460, height = 550)
dotplot(ora_analysis_cc_final, showCategory = 10)
dev.off()

#Molecular Function
ora_analysis_mf <- enrichGO(gene = ann_diff$entrezgene_id, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF",
                            pAdjustMethod = "BH", qvalueCutoff = 0.05, readable = FALSE, pool = FALSE)
ora_analysis_mf_df <- as.data.frame(ora_analysis_mf)
ora_analysis_mf_final <- clusterProfiler::simplify(ora_analysis_mf)
write_delim(x = as.data.frame(ora_analysis_mf@result), file = "Molecular Function.csv", delim = ",")
#dotplot
png("MF dotplot.png", width = 460, height = 550)
dotplot(ora_analysis_mf_final, showCategory = 10)
dev.off()

