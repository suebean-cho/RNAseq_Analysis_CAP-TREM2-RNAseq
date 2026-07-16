# ==============================================================================
# RNA-seq Analysis Pipeline for Fig 3 & Fig S5b
# ==============================================================================
library(tidyverse)
library(edgeR) 
library(amap)
library(gplots)
library(RColorBrewer)
library(ggrepel)
library(openxlsx)
library(clusterProfiler)
library(org.Hs.eg.db)   
library(enrichplot)  
library(here) 


# Initialize project directories
dir.create(here("plot"), showWarnings = FALSE)
dir.create(here("output"), showWarnings = FALSE)

data_file <- here("input", "gene_count fig3 figS5b.xls")
raw_data <- read_tsv(data_file, col_types = cols())

# Filter for protein-coding genes
raw_data_coding <- raw_data %>%
  filter(gene_biotype == "protein_coding") 

# ==============================================================================
# 1. edgeR Setup & DGEList 
# ==============================================================================
count_mat <- raw_data_coding %>%
  dplyr::select(gene_id, starts_with("Control"), starts_with("Infection")) %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix()

gene_names <- raw_data_coding %>%
  dplyr::select(gene_id, gene_name) %>%
  deframe()

gene_lengths <- raw_data_coding %>%
  dplyr::select(gene_id, gene_length) %>%
  deframe()

group <- factor(c(rep("control", 3), rep("infection", 3)), levels = c("control", "infection"))

dge <- DGEList(counts = count_mat, group = group)

# Filter low-expressed genes
keep_genes <- filterByExpr(dge)
dge <- dge[keep_genes, , keep.lib.sizes = FALSE]  
count_mat <- count_mat[keep_genes,]
gene_lengths <- gene_lengths[keep_genes]
gene_names <- gene_names[keep_genes]

# ==============================================================================
# 2. Normalization & Dispersion
# ==============================================================================
dge <- calcNormFactors(dge)
design <- model.matrix(~ 0 + group)
dge <- estimateDisp(dge, design)

plotMDS(dge)

# ==============================================================================
# 3. GLM & DEG Analysis (Infection vs Control)
# ==============================================================================
glm_fit <- glmFit(dge, design)

cont_matrix <- makeContrasts(
  Inf_vs_Ctrl = groupinfection - groupcontrol,
  levels = design
)

lrt_inf_vs_ctrl <- glmLRT(glm_fit, contrast = cont_matrix[, "Inf_vs_Ctrl"])

res_inf_vs_ctrl <- topTags(lrt_inf_vs_ctrl, n = Inf) %>% 
  as.data.frame() %>% 
  arrange(desc(logFC), FDR)

# Extract significant ENSEMBL IDs
up_genes_ensembl <- res_inf_vs_ctrl %>%
  filter(logFC > 1 & FDR < 0.05) %>%
  rownames() 

down_genes_ensembl <- res_inf_vs_ctrl %>%
  filter(logFC < -1 & FDR < 0.05) %>%
  rownames() 

# ==============================================================================
# 4. Normalized Value (TPM) Calculation
# ==============================================================================
gene_lengths_kb <- gene_lengths / 1000
rpk_mat <- dge$counts / gene_lengths_kb
tpm_mat <- t(t(rpk_mat) / colSums(rpk_mat)) * 1000000

log2_tpm_df <- as.data.frame(log2(tpm_mat + 1))

# ==============================================================================
# 5. Fig 3: Heatmap Generation
# ==============================================================================
target_genes_ensembl <- unique(c(up_genes_ensembl, down_genes_ensembl))

target_tpm_df <- log2_tpm_df %>%
  rownames_to_column(var = "gene_id") %>%          
  filter(gene_id %in% target_genes_ensembl) %>%         
  arrange(match(gene_id, target_genes_ensembl)) %>%     
  mutate(gene_name = gene_names[gene_id]) %>%          
  mutate(unique_name = make.unique(gene_name)) %>% 
  dplyr::select(-gene_id, -gene_name) %>%                  
  column_to_rownames(var = "unique_name")          

heatmap_mat <- as.matrix(target_tpm_df)
col_colors <- c(rep("#95a5a6", 3), rep("#d35400", 3))
names(col_colors) <- colnames(heatmap_mat)

heatmapCol <- function(data, col, lim, na.rm = TRUE){
  data.range <- range(data, na.rm = na.rm) 
  if(diff(data.range) == 0) stop("data has range 0")
  if(lim <= 0) stop("lim has to be positive")
  nrcol <- length(col)
  rang <- data.range[2] - data.range[1]
  reps2 <- ceiling(lim * nrcol / rang)
  col1 <- c(rep(col[1], reps2), col, rep(col[nrcol], reps2))
  return(col1)
}

rd_bu_continuous <- colorRampPalette(brewer.pal(11, "RdBu"))
heatmap_palette <- heatmapCol(data = heatmap_mat, col = rev(rd_bu_continuous(256)), lim = 2) 

hclust_func <- function(n) hcluster(n, method="pearson", link="ward", nbproc=1, doubleprecision=TRUE)

tiff(here("plot", "fig3.heatmap.tiff"), width = 16, height = 10, units = "in", res = 300, compression = "lzw") 

par(lwd = 2, cex.main = 1.5)

heatmap_obj <- heatmap.2(heatmap_mat,
                         hclustfun = hclust_func,
                         Rowv = NA,           
                         dendrogram = "col",  
                         useRaster = TRUE,    
                         symbreak = TRUE, symm = FALSE, 
                         scale = "row",       
                         density.info = "none", 
                         trace = "none", 
                         col = heatmap_palette,
                         key = TRUE, keysize = 1.0, 
                         key.par = list(cex.lab = 2.0, cex.axis = 1.5),
                         labRow = FALSE,
                         cexCol = 2.5,      
                         cex.main = 2,    
                         margin = c(12, 15), 
                         lwid = c(1, 4),  
                         lhei = c(1, 5),  
                         main = "\n",
                         ColSideColors = col_colors)

par(xpd=TRUE)
legend("topright", 
       legend = c("Control", "Infection"), 
       fill = c("#95a5a6", "#d35400"),     
       bty = "n",                    
       cex = 2.0,                    
       inset = c(0, -0.05))          

dev.off()

# ==============================================================================
# 6. Save Heatmap Source Data (Supplementary Table)
# ==============================================================================
heatmap_stats_df <- res_inf_vs_ctrl %>%
  rownames_to_column(var = "gene_id") %>%
  filter(gene_id %in% target_genes_ensembl) %>%
  mutate(gene_name = gene_names[gene_id]) %>%
  mutate(Direction = ifelse(logFC > 0, "UP", "DOWN")) %>% 
  dplyr::select(gene_id, gene_name, Direction, logFC, logCPM, LR, PValue, FDR)

raw_counts_df <- as.data.frame(count_mat[target_genes_ensembl, ])
colnames(raw_counts_df) <- paste0("Raw_", colnames(raw_counts_df))
raw_counts_df <- raw_counts_df %>% rownames_to_column(var = "gene_id")

tpm_export_df <- as.data.frame(tpm_mat[target_genes_ensembl, ])
colnames(tpm_export_df) <- paste0("TPM_", colnames(tpm_export_df))
tpm_export_df <- tpm_export_df %>% rownames_to_column(var = "gene_id")

heatmap_supplementary_df <- heatmap_stats_df %>%
  left_join(raw_counts_df, by = "gene_id") %>%
  left_join(tpm_export_df, by = "gene_id") %>%
  arrange(desc(logFC)) 

write.xlsx(heatmap_supplementary_df, 
           file = here("output", "fig3.Heatmap_Source_Data.xlsx"), 
           rowNames = FALSE)

# ==============================================================================
# 7. Fig 3: Volcano Plot
# ==============================================================================
volcano_df <- res_inf_vs_ctrl %>%
  rownames_to_column(var = "gene_id") %>%          
  mutate(gene_name = gene_names[gene_id]) 

volcano_df$diffexpressed <- "NO"
volcano_df$diffexpressed[volcano_df$logFC > 1 & volcano_df$FDR < 0.05] <- "UP"
volcano_df$diffexpressed[volcano_df$logFC < -1 & volcano_df$FDR < 0.05] <- "DOWN"

target_genes_volcano <- volcano_df %>% 
  filter(gene_name %in% c("CCL4","CCL2","CISH","IL23A","EBI3","PTGS2", "PTGES", "IL6", "IL1B"))

tiff(here("plot", "fig3.volcanoplot.tiff"), width = 8, height = 8, units = "in", res = 300, compression = "lzw") 

plot_volcano <- ggplot(data = volcano_df) + aes(x=logFC, y=-log10(FDR)) +
  geom_point(aes(colour = diffexpressed), alpha = 0.7, shape = 16, size = 2.5, na.rm = TRUE) + 
  theme_minimal() +
  geom_text_repel(data = target_genes_volcano, 
                  aes(label = gene_name),
                  force = 2,
                  nudge_y = 10,
                  size = 7,
                  fontface = "bold", 
                  min.segment.length = 0,      
                  segment.color = "black",     
                  segment.size = 0.5) +        
  geom_vline(xintercept = c(-1,1), col = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = 'black', linetype = 'dashed') +
  coord_cartesian(xlim = c(-10,10), ylim = c(0,190)) +
  scale_color_manual(values = c('royalblue','grey','red2')) +
  theme(text = element_text(size = 22, face = 'bold'),     
        axis.title = element_text(size = 24),   
        axis.text.x = element_text(size = 18, color = "black"), 
        axis.text.y = element_text(size = 18, color = "black"),
        legend.position = "none",
        axis.line.x.bottom = element_line(linetype="solid", size=1),
        axis.line.y.left = element_line(linetype="solid", size=1)) +
  labs(x = expression("Log"["2"]*" FC"), y = expression("-Log"["10"]*" FDR"))

print(plot_volcano)
dev.off()

# ==============================================================================
# 8. Fig S5b: ORA (KEGG & GO-BP) Setup
# ==============================================================================
bg_ensembl <- rownames(res_inf_vs_ctrl)

bg_map <- bitr(bg_ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) 
bg_entrez <- unique(bg_map$ENTREZID) 

up_map <- bitr(up_genes_ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) 
up_entrez <- unique(up_map$ENTREZID)

down_map <- bitr(down_genes_ensembl, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) 
down_entrez <- unique(down_map$ENTREZID) 

# ==============================================================================
# 9. Function: Run ORA Pipeline (KEGG & GO-BP)
# ==============================================================================
run_ora_pipeline <- function(entrez_genes, bg_genes, db_type, x_limit, plot_title, file_name) {
  
  # Enrichment Analysis
  if (db_type == "KEGG") {
    ora_obj <- enrichKEGG(gene = entrez_genes, organism = "hsa", pvalueCutoff = 1.0, qvalueCutoff = 1.0, universe = bg_genes)
  } else {
    ont_type <- sub("GO_", "", db_type)
    ora_obj <- enrichGO(gene = entrez_genes, OrgDb = org.Hs.eg.db, ont = ont_type, pvalueCutoff = 1.0, qvalueCutoff = 1.0, universe = bg_genes)
  }
  # Extract top 25 for plotting
  ora_top25 <- ora_obj
  ora_top25@result <- ora_top25@result %>% 
    arrange(p.adjust) %>% 
    slice_head(n = 25)
  
  # Plot dotplot (p.adjust > 0.05 will map to NA/grey via scale_fill_gradient)
  p <- dotplot(ora_top25, showCategory = 25, orderBy = "p.adjust", decreasing = FALSE) +
    scale_fill_gradient(limits = c(0, 0.05), low = "red", high = "blue", na.value = "grey50") +
    scale_x_continuous(limits = c(0, x_limit)) + 
    ggtitle(plot_title) +
    theme(text = element_text(size = 18),
          axis.text.y = element_text(size = 14),
          axis.text.x = element_text(size = 16),
          plot.title = element_text(size = 18, face = "bold", hjust = 0.5)) +
    guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2), fill = guide_colorbar(order = 2))
  
  ggsave(filename = here("plot", file_name), plot = p, width = 8, height = 12, units = "in", dpi = 300, compression = "lzw")
  
  # Export configuration: all significant terms (p < 0.05) if > 25, otherwise top 25
  ora_readable <- setReadable(ora_obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  all_res_df <- ora_readable@result %>% arrange(p.adjust)
  sig_count <- sum(all_res_df$p.adjust < 0.05, na.rm = TRUE)
  
  if (sig_count > 25) {
    res_df <- all_res_df %>% filter(p.adjust < 0.05)
  } else {
    res_df <- all_res_df %>% slice_head(n = 25)
  }
  
  return(res_df)
}

# ==============================================================================
# 10. Execute ORA Pipeline and Export Results
# ==============================================================================
# KEGG Analysis
kegg_up_res <- run_ora_pipeline(up_entrez, bg_entrez, "KEGG", 0.17, 
                                "Up-regulated KEGG Pathway enrichment (Top 25)", "figS5b.kegg_up_dotplot.tiff")
kegg_down_res <- run_ora_pipeline(down_entrez, bg_entrez, "KEGG", 0.10, 
                                  "Down-regulated KEGG Pathway enrichment (Top 25)", "figS5b.kegg_down_dotplot.tiff")

write.xlsx(list("Up-regulated" = kegg_up_res, "Down-regulated" = kegg_down_res), 
           file = here("output", "figS5b.KEGG_Pathway_Results.xlsx"), rowNames = FALSE)

# GO_BP Analysis
go_bp_up_res <- run_ora_pipeline(up_entrez, bg_entrez, "GO_BP", 0.15, 
                                 "Up-regulated GO Biological Process (Top 25)", "figS5b.go_bp_up_dotplot.tiff")
go_bp_down_res <- run_ora_pipeline(down_entrez, bg_entrez, "GO_BP", 0.10, 
                                   "Down-regulated GO Biological Process (Top 25)", "figS5b.go_bp_down_dotplot.tiff")

write.xlsx(list("Up-regulated" = go_bp_up_res, "Down-regulated" = go_bp_down_res), 
           file = here("output", "figS5b.GO_BP_Pathway_Results.xlsx"), rowNames = FALSE)

# GO_CC Analysis
go_cc_up_res <- run_ora_pipeline(up_entrez, bg_entrez, "GO_CC", 0.15, 
                                 "Up-regulated GO Cellular Component (Top 25)", "figS5b.go_cc_up_dotplot.tiff")
go_cc_down_res <- run_ora_pipeline(down_entrez, bg_entrez, "GO_CC", 0.10, 
                                   "Down-regulated GO Cellular Component (Top 25)", "figS5b.go_cc_down_dotplot.tiff")

write.xlsx(list("Up-regulated" = go_cc_up_res, "Down-regulated" = go_cc_down_res), 
           file = here("output", "figS5b.GO_CC_Pathway_Results.xlsx"), rowNames = FALSE)

# GO_MF Analysis
go_mf_up_res <- run_ora_pipeline(up_entrez, bg_entrez, "GO_MF", 0.15, 
                                 "Up-regulated GO Molecular Function (Top 25)", "figS5b.go_mf_up_dotplot.tiff")
go_mf_down_res <- run_ora_pipeline(down_entrez, bg_entrez, "GO_MF", 0.10, 
                                   "Down-regulated GO Molecular Function (Top 25)", "figS5b.go_mf_down_dotplot.tiff")

write.xlsx(list("Up-regulated" = go_mf_up_res, "Down-regulated" = go_mf_down_res), 
           file = here("output", "figS5b.GO_MF_Pathway_Results.xlsx"), rowNames = FALSE)

# ==============================================================================
# 11. Reproducibility Logging
# ==============================================================================
writeLines(capture.output(sessionInfo()), here("output", "sessionInfo.txt"))