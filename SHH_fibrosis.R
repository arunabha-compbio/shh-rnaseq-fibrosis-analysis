# ----------------------------------------------------------------------------- 
# STEP 1: Install Required Packages (Run ONCE, then comment out) 
# ----------------------------------------------------------------------------- 
install.packages(c("ggplot2", "pheatmap", "ggrepel")) 
if (!require("BiocManager")) install.packages("BiocManager") 
BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot")) 
# ----------------------------------------------------------------------------- 
# STEP 2: Load Libraries 
# ----------------------------------------------------------------------------- 
library(ggplot2) 
library(pheatmap) 
library(ggrepel) 
library(clusterProfiler) 
library(org.Mm.eg.db) 
library(enrichplot) 
# ----------------------------------------------------------------------------- 
# STEP 3: Load the Count Data 
# 
# Your file is a featureCounts output. 
# - First row  : comment line (starts with #)  → skipped automatically 
# - Second row : column headers 
# - Columns 1-6: Gene info (Geneid, Chr, Start, End, Strand, Length) 
# - Columns 7+ : Sample counts (Control, GR1) 
# ----------------------------------------------------------------------------- 
counts <- read.table("counts.txt", 
                     header        
                     = TRUE, 
                     sep           = "\t", 
                     comment.char  = "#",   # skips the first comment line 
                     check.names   = FALSE) 
# Set gene IDs as row names 
rownames(counts) <- counts$Geneid 
# Extract ONLY the count columns (columns 7 onward) 
counts_clean <- counts[, 7:ncol(counts)] 
# Rename columns to simple sample names 
colnames(counts_clean) <- c("Control", "GR1") 
# Quick check 
head(counts_clean) 
dim(counts_clean)   # shows total number of genes 
# ----------------------------------------------------------------------------- 
# STEP 4: Filter Low-Count Genes 
# Remove genes where total counts across all samples is <= 10 
# (These are likely noise or unexpressed genes) 
# ----------------------------------------------------------------------------- 
counts_clean <- counts_clean[rowSums(counts_clean) > 10, ] 
cat("Genes remaining after filtering:", nrow(counts_clean), "\n") 
# ----------------------------------------------------------------------------- 
# STEP 5: Log2 Transformation 
# 
# We add 1 before log to avoid log(0) = -Inf 
# This is called "log2(counts + 1)" or "log2CPM-like" transformation 
# ----------------------------------------------------------------------------- 
log_counts <- log2(counts_clean + 1) 
# ----------------------------------------------------------------------------- 
# STEP 6: Calculate Fold Change (GR1 vs Control) 
# 
# logFC = log2(GR1) - log2(Control) 
# Positive logFC = higher in GR1 (upregulated) 
# Negative logFC = lower in GR1 (downregulated) 
# ----------------------------------------------------------------------------- 
logFC <- log_counts[, "GR1"] - log_counts[, "Control"] 
res <- data.frame( 
  ENSEMBL        = rownames(log_counts), 
  log2FoldChange = logFC 
) 
# Sort by fold change (highest first) 
res <- res[order(res$log2FoldChange, decreasing = TRUE), ] 
head(res) 
# ----------------------------------------------------------------------------- 
# STEP 7: Convert ENSEMBL IDs → Gene Symbols 
# 
# Your gene IDs look like: ENSMUSG00000029844.16 
# We must remove the version number (.16) before mapping 
# ----------------------------------------------------------------------------- 
res$ENSEMBL_clean <- sub("[.].*", "", res$ENSEMBL) 
gene_map <- bitr(res$ENSEMBL_clean, 
                 fromType = "ENSEMBL", 
                 toType   = c("SYMBOL", "ENTREZID"), 
                 OrgDb    = org.Mm.eg.db) 
# Merge symbol information into results 
res_final <- merge(res, gene_map, 
                   by.x = "ENSEMBL_clean", 
                   by.y = "ENSEMBL") 
# Clean up 
res_final <- res_final[!is.na(res_final$SYMBOL), ] 
res_final <- res_final[res_final$SYMBOL != "", ] 
res_final <- res_final[!duplicated(res_final$SYMBOL), ] 
# Keep useful columns 
res_final <- res_final[, c("SYMBOL", "log2FoldChange", "ENTREZID", 
                           "ENSEMBL_clean")] 
rownames(res_final) <- res_final$SYMBOL 
cat("Total genes with symbol:", nrow(res_final), "\n") 
head(res_final) 
# 
==================================================================
  =========== 
  #  SHH PATHWAY ANALYSIS 
  # 
  #  Sonic Hedgehog (SHH) pathway components in mouse: 
  # 
  #  LIGANDS    : Shh, Ihh, Dhh 
  #  RECEPTORS  : Ptch1, Ptch2 (inhibitory), Smo (activating) 
  #  CO-RECEPT. : Gas1, Cdon, Boc 
  #  SIGNAL     : Sufu (repressor), Kif7 
  #  EFFECTORS  : Gli1, Gli2, Gli3 (transcription factors) 
  #  FEEDBACK   : Hhip (inhibitor), Ptch1 (feedback repressor) 
  #  TARGETS    : Ccnd1, Mycn, Snai1, Foxf1 
  # 
  ==================================================================
  =========== 
  # ----------------------------------------------------------------------------- 
# STEP 8: Define SHH Pathway Gene List 
# ----------------------------------------------------------------------------- 
shh_genes <- c( 
  # Ligands 
  "Shh", "Ihh", "Dhh", 
  # Receptors & Co-receptors 
  "Ptch1", "Ptch2", "Smo", 
  "Gas1", "Cdon", "Boc", 
  # Intracellular signal transducers 
  "Sufu", "Kif7", 
  # GLI transcription factors (main effectors) 
  "Gli1", "Gli2", "Gli3", 
  # Feedback regulators & targets 
  "Hhip", 
  "Ccnd1",   # Cyclin D1 - cell cycle target 
  "Mycn",    # N-Myc - proliferation target 
  "Snai1",   # Snail - EMT / fibrosis 
  "Foxf1"    # Forkhead - mesenchymal marker 
) 
# ----------------------------------------------------------------------------- 
# STEP 9: Extract SHH Genes from Results 
# ----------------------------------------------------------------------------- 
shh_data <- res_final[res_final$SYMBOL %in% shh_genes, ] 
shh_data <- shh_data[order(shh_data$log2FoldChange, decreasing = TRUE), ] 
cat("\n--- SHH Pathway Genes Found in Dataset ---\n") 
print(shh_data[, c("SYMBOL", "log2FoldChange")]) 
# How many were detected? 
cat("\nTotal SHH genes detected:", nrow(shh_data), "out of", length(shh_genes), "\n") 
# ----------------------------------------------------------------------------- 
# STEP 10: Determine Activation or Suppression 
# 
# logFC > 0.5  → Upregulated / ACTIVATED in GR1 
# logFC < -0.5 → Downregulated / SUPPRESSED in GR1 
# -0.5 to 0.5  → No significant change 
# ----------------------------------------------------------------------------- 
shh_data$Status <- ifelse(shh_data$log2FoldChange >  0.5, "Activated", 
                          ifelse(shh_data$log2FoldChange < -0.5, "Suppressed", 
                                 "No Change")) 
cat("\n--- Activation/Suppression Status ---\n") 
print(shh_data[, c("SYMBOL", "log2FoldChange", "Status")]) 
# ----------------------------------------------------------------------------- 
# STEP 11: Bar Plot – SHH Gene Expression (GR1 vs Control) 
# ----------------------------------------------------------------------------- 
ggplot(shh_data, 
       aes(x    = reorder(SYMBOL, log2FoldChange), 
           y    = log2FoldChange, 
           fill = log2FoldChange)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  coord_flip() + 
  scale_fill_gradient2(low
                       =
                         "steelblue", 
                       mid      = "white", 
                       high     = "firebrick", 
                       midpoint = 0, 
                       name     = "log2FC") + 
  geom_hline(yintercept =  0.5, linetype = "dashed", color = "gray40") + 
  geom_hline(yintercept = -0.5, linetype = "dashed", color = "gray40") + 
  theme_minimal(base_size = 13) + 
  labs(title    = "SHH Signaling Pathway – Gene Expression", 
       subtitle = "GR1 vs Control | Red = Upregulated, Blue = Downregulated", 
       x        = "Gene", 
       y        = "log2 Fold Change") + 
  theme(plot.title    = element_text(face = "bold"), 
        plot.subtitle = element_text(color = "gray50")) 
# ----------------------------------------------------------------------------- 
# STEP 12: Heatmap – SHH Gene Expression Pattern 
# ----------------------------------------------------------------------------- 
# Get log-counts for SHH genes (both samples) 
# First, merge log_counts with gene symbols 
log_counts_df <- as.data.frame(log_counts) 
log_counts_df$ENSEMBL_clean <- sub("\\..*", "", rownames(log_counts_df)) 
log_merged <- merge(log_counts_df, gene_map, 
                    by.x = "ENSEMBL_clean", 
                    by.y = "ENSEMBL") 
log_merged <- log_merged[!is.na(log_merged$SYMBOL), ] 
log_merged <- log_merged[!duplicated(log_merged$SYMBOL), ] 
rownames(log_merged) <- log_merged$SYMBOL 
log_final <- log_merged[, c("Control", "GR1")] 
# Filter for SHH genes present in data 
shh_present <- shh_genes[shh_genes %in% rownames(log_final)] 
shh_matrix  <- log_final[shh_present, ] 
cat("\nSHH genes in heatmap:", nrow(shh_matrix), "\n") 
# Add annotation for activation status 
annotation_row <- data.frame( 
  Status = ifelse(shh_data[shh_present, "log2FoldChange"] >  0.5, "Activated", 
                  ifelse(shh_data[shh_present, "log2FoldChange"] < -0.5, "Suppressed", 
                         "No Change")) 
) 
rownames(annotation_row) <- shh_present 
ann_colors <- list( 
  Status = c(Activated  = "firebrick", 
             Suppressed = "steelblue", 
             `No Change` = "gray80") 
) 
pheatmap(shh_matrix, 
         scale          = "row",       
         cluster_rows   = TRUE, 
         cluster_cols   = FALSE, 
         # z-score per gene 
         annotation_row = annotation_row, 
         annotation_colors = ann_colors, 
         color          = colorRampPalette(c("steelblue","white","firebrick"))(100), 
         fontsize_row   = 11, 
         fontsize_col   = 12, 
         main           = "SHH Pathway Gene Expression Heatmap\n(Row-scaled, GR1 vs 
Control)") 
# ----------------------------------------------------------------------------- 
# STEP 13: Bubble / Dot Plot – Effect Size Summary 
# 
# This shows each SHH gene's expression level (x-axis) and 
# fold-change magnitude (bubble size) 
# ----------------------------------------------------------------------------- 
shh_plot_df <- shh_data 
shh_plot_df$MeanExpr <- rowMeans(shh_matrix[shh_plot_df$SYMBOL, ], na.rm = TRUE) 
shh_plot_df$AbsFC    <- abs(shh_plot_df$log2FoldChange) 
ggplot(shh_plot_df, 
       aes(x     = MeanExpr, 
           y     = log2FoldChange, 
           size  = AbsFC, 
           color = Status, 
           label = SYMBOL)) + 
  geom_point(alpha = 0.8) + 
  geom_text_repel(size = 3.5, max.overlaps = 15) + 
  scale_color_manual(values = c("Activated"  = "firebrick", 
                                "Suppressed" = "steelblue", 
                                "No Change"  = "gray60")) + 
  scale_size_continuous(range = c(3, 10)) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") + 
  theme_minimal(base_size = 13) + 
  labs(title    = "SHH Pathway – Expression vs Fold Change", 
       subtitle = "Bubble size = magnitude of change", 
       x        = "Mean log2 Expression", 
       y        = "log2 Fold Change (GR1 vs Control)", 
       color    = "Status", 
       size
       = "|log2FC|") 
# ----------------------------------------------------------------------------- 
# STEP 14: Pathway Role Summary Plot (Grouped by function) 
# ----------------------------------------------------------------------------- 
shh_data$Category <- ifelse(shh_data$SYMBOL %in% c("Shh","Ihh","Dhh"), 
                            "Ligand", 
                            ifelse(shh_data$SYMBOL %in% c("Ptch1","Ptch2","Smo", 
                                                          "Gas1","Cdon","Boc"), 
                                   "Receptor/Co-receptor", 
                                   ifelse(shh_data$SYMBOL %in% c("Gli1","Gli2","Gli3"), 
                                          "GLI Transcription Factor", 
                                          ifelse(shh_data$SYMBOL %in% c("Sufu","Kif7","Hhip"), 
                                                 "Regulator", 
                                                 "Downstream Target")))) 
ggplot(shh_data, 
       aes(x    = reorder(SYMBOL, log2FoldChange), 
           y    = log2FoldChange, 
           fill = Category)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  coord_flip() + 
  scale_fill_manual(values = c( 
    "Ligand"                    = "#E67E22", 
    "Receptor/Co-receptor"      = "#8E44AD", 
    "GLI Transcription Factor"  = "#E74C3C", 
    "Regulator"                 = "#27AE60", 
    "Downstream Target"  
  )) + 
  =
  "#2980B9" 
geom_hline(yintercept = 0, color = "black", linewidth = 0.5) + 
  theme_minimal(base_size = 13) + 
  labs(title    = "SHH Pathway – Gene Expression by Functional Category", 
       subtitle = "Positive = upregulated in GR1 | Negative = downregulated in GR1", 
       x        = "Gene", 
       y        = "log2 Fold Change", 
       fill    
       = "Pathway Role") 
# ----------------------------------------------------------------------------- 
# STEP 15: Fibrosis & Repair Connection 
# 
# SHH drives fibrosis through GLI2 → TGF-β / CTGF axis 
# We check expression of key fibrosis genes alongside SHH effectors 
# ----------------------------------------------------------------------------- 
fibrosis_repair_genes <- c( 
  # SHH effectors (fibrosis-relevant) 
  "Gli1", "Gli2", "Smo", "Ptch1", 
  # TGF-β / fibrosis axis 
  "Tgfb1", "Tgfb2", "Tgfbr1", 
  # Extracellular matrix / scar genes 
  "Col1a1", "Col1a2", "Col3a1", "Fn1", 
  # Myofibroblast marker 
  "Acta2", 
  # Repair / regeneration 
  "Vegfa", "Mki67", "Cdh1" 
) 
# Extract from full results 
fibro_data <- res_final[res_final$SYMBOL %in% fibrosis_repair_genes, ] 
fibro_data <- fibro_data[order(fibro_data$log2FoldChange, decreasing = TRUE), ] 
cat("\n--- Fibrosis & Repair Gene Expression ---\n") 
print(fibro_data[, c("SYMBOL", "log2FoldChange")]) 
# Plot 
ggplot(fibro_data, 
       aes(x    = reorder(SYMBOL, log2FoldChange), 
           y    = log2FoldChange, 
           fill = log2FoldChange)) + 
  geom_bar(stat = "identity", width = 0.7) + 
  coord_flip() + 
  scale_fill_gradient2(low
                       =
                         "steelblue", 
                       mid      = "white", 
                       high     = "firebrick", 
                       midpoint = 0) + 
  geom_hline(yintercept = 0, color = "black") + 
  theme_minimal(base_size = 13) + 
  labs(title    = "SHH Pathway Effectors & Fibrosis/Repair Genes", 
       subtitle = "Connecting SHH signaling to tissue remodeling", 
       x        = "Gene", 
       y        = "log2 Fold Change (GR1 vs Control)", 
       fill    
       = "log2FC") 
# ----------------------------------------------------------------------------- 
# STEP 16: GSEA – GO Biological Process Enrichment (full dataset) 
# ----------------------------------------------------------------------------- 
# Build ranked gene list using all expressed genes 
geneList <- res_final$log2FoldChange 
names(geneList) <- res_final$ENTREZID 
geneList <- geneList[!is.na(names(geneList))] 
geneList <- sort(geneList, decreasing = TRUE) 
# Run GSEA on GO Biological Process 
gsea_go <- gseGO(geneList = geneList, 
                 OrgDb    = org.Mm.eg.db, 
                 ont      = "BP", 
                 minGSSize = 10, 
                 pvalueCutoff = 0.05) 
# View top enriched pathways 
head(as.data.frame(gsea_go)[, c("Description", "NES", "pvalue")], 15) 
# Dot plot of top enriched GO terms 
dotplot(gsea_go, showCategory = 15, 
        title = "GSEA – GO Biological Process Enrichment") 
# ----------------------------------------------------------------------------- 
# STEP 17: KEGG Pathway Enrichment 
# ----------------------------------------------------------------------------- 
gsea_kegg <- gseKEGG(geneList = geneList, 
                     organism = "mmu",         # mmu = Mus musculus (mouse) 
                     minGSSize = 10, 
                     pvalueCutoff = 0.05) 
dotplot(gsea_kegg, showCategory = 10, 
        title = "GSEA – KEGG Pathway Enrichment (Mouse)") 


