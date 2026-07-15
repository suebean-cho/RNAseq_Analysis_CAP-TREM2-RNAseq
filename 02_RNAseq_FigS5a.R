# ==============================================================================
# RNA-seq Analysis Pipeline for Fig S5a
# ==============================================================================
library(tidyverse)
library(edgeR) 
library(amap)
library(gplots)
library(RColorBrewer)
library(ggrepel)
library(openxlsx)
library(pheatmap)  
library(here) 

# Initialize project directories
dir.create(here("plot"), showWarnings = FALSE)
dir.create(here("output"), showWarnings = FALSE)

data_file <- here("input", "gene_count figS5a.xls")
raw_data <- read_tsv(data_file, col_types = cols())

# Filter for protein-coding genes
raw_data_coding <- raw_data %>%
  filter(gene_biotype == "protein_coding") 

# ==============================================================================
# 1. edgeR Setup & DGEList 
# ==============================================================================
count_mat <- raw_data_coding %>%
  dplyr::select(gene_id, Control, CA_D1, CA_D2, CA_C1) %>%
  column_to_rownames(var = "gene_id") %>%
  as.matrix()

gene_names <- raw_data_coding %>%
  dplyr::select(gene_id, gene_name) %>%
  deframe()

gene_lengths <- raw_data_coding %>%
  dplyr::select(gene_id, gene_length) %>%
  deframe()

group <- factor(c("Control", rep("CA_G2", 3)), levels = c("Control", "CA_G2"))

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
design <- model.matrix(~ group)
dge <- estimateDisp(dge, design)

col_sub <- c("blue", rep("red",3))
plotMDS(dge, col = col_sub)

# ==============================================================================
# 3. GLM & DEG Analysis (CA_G2 vs Control)
# ==============================================================================
glm_fit <- glmFit(dge, design)

lrt_inf_vs_ctrl <- glmLRT(glm_fit, contrast = c(0, 1))

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
# 5. Fig S5a: Heatmap Generation
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
col_colors <- c(rep("#95a5a6", 1), rep("#d35400", 3))
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

hclust_func <- function(n) hcluster(n, method="pearson", link="ward", nbproc=4, doubleprecision=TRUE)

tiff(here("plot", "figS5a.heatmap.tiff"), width = 16, height = 10, units = "in", res = 300, compression = "lzw") 

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
                         adjCol = c(1, 0.5),
                         cex.main = 2,    
                         margin = c(12, 15), 
                         lwid = c(1, 4),  
                         lhei = c(1, 5),  
                         main = "\n",
                         ColSideColors = col_colors)

par(xpd=TRUE)
legend("topright", 
       legend = c("Control", "CA_G2"), 
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
           file = here("output", "figS5a.Heatmap_Source_Data.xlsx"), 
           rowNames = FALSE)

# ==============================================================================
# 7. Fig S5a: Volcano Plot
# ==============================================================================
volcano_df <- res_inf_vs_ctrl %>%
  rownames_to_column(var = "gene_id") %>%          
  mutate(gene_name = gene_names[gene_id]) 

volcano_df$diffexpressed <- "NO"
volcano_df$diffexpressed[volcano_df$logFC > 1 & volcano_df$FDR < 0.05] <- "UP"
volcano_df$diffexpressed[volcano_df$logFC < -1 & volcano_df$FDR < 0.05] <- "DOWN"

target_genes_volcano <- volcano_df %>% 
  filter(gene_name %in% c("Il6", "Il1b", "Il1a", "Ccl7", "Ccl2"))

tiff(here("plot", "figS5a.volcanoplot.tiff"), width = 8, height = 8, units = "in", res = 300, compression = "lzw") 

plot_volcano <- ggplot(data = volcano_df) + aes(x=logFC, y=-log10(FDR)) +
  geom_point(aes(colour = diffexpressed), alpha = 0.7, shape = 16, size = 2.5, na.rm = TRUE) + 
  theme_minimal() +
  geom_text_repel(data = target_genes_volcano, 
                  aes(label = gene_name),
                  force = 2,
                  nudge_y = 5,
                  size = 7,
                  fontface = "bold", 
                  min.segment.length = 0,      
                  segment.color = "black",     
                  segment.size = 0.5) +        
  geom_vline(xintercept = c(-1,1), col = 'black', linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = 'black', linetype = 'dashed') +
  coord_cartesian(xlim = c(-16,16), ylim = c(0,150)) +
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
# 8. Fig S5a: Pairs Plot (Scatter)
# ==============================================================================
data_g2 <- log2_tpm_df[, -1]

panel_cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y)
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.7/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * r)
}

tiff(filename = here("plot", "figS5a.ScatterPlot.tiff"), width = 10, height = 10, units = "in", res = 300, compression = "lzw")

pairs(data_g2, 
      lower.panel = panel.smooth, 
      upper.panel = panel_cor,     
      pch = 19, col = rgb(0,0,0,0.1), 
      main = "",
      cex.labels = 7.0,  
      cex.axis = 2.5)    

dev.off() 

# ==============================================================================
# 9. Fig S5a: Correlation Matrix
# ==============================================================================
sample_cor <- cor(log2_tpm_df) 

pheatmap(sample_cor, 
         main = "",
         filename = here("plot", "figS5a.CorrMatrix.tiff"), 
         width = 8, height = 7,               
         res = 300,                           
         fontsize = 16,          
         fontsize_row = 18,      
         fontsize_col = 18,      
         display_numbers = FALSE,       
         number_color = "black",       
         fontsize_number = 14)         

# ==============================================================================
# 10. Reproducibility Logging
# ==============================================================================
#writeLines(capture.output(sessionInfo()), here("output", "sessionInfo.txt"))