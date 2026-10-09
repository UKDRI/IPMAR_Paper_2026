# limma analysis with donor blocking and keratin exclusion. 
suppressPackageStartupMessages({
  library(readr)
  library(readxl)
  library(dplyr)
  library(limma)
  library(ggplot2)
  library(EnhancedVolcano)
  library(org.Hs.eg.db)
  library(fgsea)
})

# 0. Paths & parameters
base_dir <- "/scratch/SCWF00010/rebecca/proteomics"
outdir   <- file.path(base_dir, "analysis_results_unimputed_all_no_keratins_jun26_excluded_samples_add54A/")
dir.create(outdir, showWarnings = FALSE)

min_detect <- 2        # >=2 samples per condition
log_transform <- TRUE # assume raw intensities

# 1. Metadata
meta <- read_xlsx(file.path(base_dir, "IPMAR_Proteomics_Masterfile.xlsx"))
meta$cnames <- gsub(".*IMPAR", "IMPAR", meta$Proteomics_Internal.code)
meta <- as.data.frame(meta)
rownames(meta) <- meta$cnames

meta <- meta %>%
  mutate(
    AD_Status = factor(AD_Status, levels = c("Control","LOAD")),
    condition = factor(Proteomics_Condition, levels = c("ctrl","lps")),
    donor     = factor(Cell_line_ID_Short)
  )

# 2. DIA-NN protein matrix
pg <- read_delim(file.path(base_dir, "allplates/report.pg_matrix.tsv"), delim = "\t")

meta_cols <- c("Protein.Group", "Protein.Names", "Genes",
               "First.Protein.Description", "N.Sequences", "N.Proteotypic.Sequences")

sample_cols <- setdiff(colnames(pg), meta_cols)
clean_names <- gsub(".*IMPAR", "IMPAR", sample_cols)
colnames(pg)[match(sample_cols, colnames(pg))] <- clean_names

expr <- pg[, clean_names]
expr <- as.matrix(sapply(expr, as.numeric))
rownames(expr) <- pg$Protein.Group

# keep samples present in metadata
expr <- expr[, colnames(expr) %in% rownames(meta)]
meta <- meta[colnames(expr), ]
stopifnot(all(colnames(expr) == rownames(meta)))
vars_needed <- c("AD_Status", "condition", "donor")
keep <- complete.cases(meta[, vars_needed])

cat("Dropping", sum(!keep), "samples due to missing metadata\n")

meta   <- meta[keep, ]
expr   <- expr[, keep]

stopifnot(
  ncol(expr) == nrow(meta),
  all(colnames(expr) == rownames(meta))
)

# 3. Log2 transform
if(log_transform){
  expr <- log2(expr)
}

#normalise
expr <- sweep(expr, 2, apply(expr, 2, median, na.rm = TRUE))

# Average technical replicates by donor and condition
meta$BioID <- paste0(meta$donor, "_", meta$condition)

expr_averaged_T <- avereps(t(expr), ID = meta$BioID)
expr_averaged <- t(expr_averaged_T)

meta_averaged <- meta[!duplicated(meta$BioID), ]

# avereps reorders the columns alphabetically by BioID
expr_averaged <- expr_averaged[, meta_averaged$BioID]

stopifnot(ncol(expr_averaged) == nrow(meta_averaged))
stopifnot(all(colnames(expr_averaged) == meta_averaged$BioID))

cat("Original columns:", ncol(expr), "\n")
cat("Averaged columns:", ncol(expr_averaged), "\n")

expr <- expr_averaged
meta <- meta_averaged
# 4. Filtering (unimputed)
ctrl_ix <- meta$condition == "ctrl"
lps_ix  <- meta$condition == "lps"

keep_ctrl <- rowSums(!is.na(expr[, ctrl_ix, drop=FALSE])) >= min_detect
keep_lps  <- rowSums(!is.na(expr[, lps_ix,  drop=FALSE])) >= min_detect
keep_any  <- keep_ctrl | keep_lps

expr_f <- expr[keep_any, ]
cat("Proteins kept:", nrow(expr_f), "\n")

# Exclude repeat batch and outlier samples
keep <- meta_averaged$Cardiff_Batch_ID != "batch_repeats"
expr2 <- expr_f[, keep]
meta2 <- meta_averaged[keep, ]

# Exclude 21A control outlier
expr3 <- expr2[,-which(colnames(expr2)=="21A_ctrl")]
meta3 <- meta2[-which(meta2$BioID=="21A_ctrl"),]
# Exclude donor 40A
expr3 <- expr3[,-grep("40A",colnames(expr3))]
meta3 <- meta3[-grep("40A",meta3$donor),]

# Strip everything after the first semicolon
rownames(expr3) <- sub(";.*", "", rownames(expr3))
expr3 <- avereps(expr3) 

# Keratin and KRTAP UniProt exclusions
keratin_ids <- c(
  # The Troublesome Unreviewed / TrEMBL Variants
  "Q5XKE5", "A0A024R254", "A0A024R314", "A0A024R1R7", "A0A024R1V0", "A0A0B4J2D5",
  "E9PH74", "F8VP11", "H0Y2A9", "G3V3H5", "V9GZA4", "Q53XF5", "A0A140VJE5",
  "A0A024R267", "A0A024R1X5", "A0A024R204", "A0A087X1V0", "A0A024R213",
  
  # Canonical Type I Cytoskeletal / Epithelial Keratins
  "P13645", "P02533", "P08779", "P48668", "P05787", "P08727", "P35900", 
  "P13646", "P19013", "Q7Z3Z6", "Q7Z3Z7", "Q7Z3Y9", "Q7Z3Z5", "Q7Z3Z4", 
  "Q7Z3Z3", "A6NCE7", "A6NCI5", "Q14525",
  
  # Canonical Type II Cytoskeletal / Epithelial Keratins
  "P04264", "P35908", "P13647", "P19012", "P35527", "P05783", "Q5XPE7", 
  "O76011", "Q6T528", "Q9C075", "P12035", "P02538", "A1X935", "Q7Z2K1", 
  "P48669", "Q14194", "Q15323", "Q14204", "Q99455", "Q8N1N4",
  
  # Hair & Cuticular Keratins (Type I & II)
  "Q92764", "O76009", "O76013", "Q14523", "Q14532", "Q6A162", "Q6A163", 
  "A8MWD9", "A1A5D5", "A8MTI9", "Q01546", "Q14533", "Q99456", "P78385", 
  "O76015", "P78386", "O43790", "Q9BYP7", "Q92765", "Q3SY84", "Q14554", 
  "Q3SY78", "Q3SXZ3", "A4D2H0", "A6NKV7",
  
  # Inner Root Sheath Keratins
  "Q96I24", "A6NLC5", "A6NLC2", "Q86Y46", "A6NFT2", "Q8IWU8", "A6NMD2", 
  "Q8N6L9", "A8MX37", "A6NJS3", "A6NH52", "A6NJT0",
  
  # Common Keratin-Associated Proteins (KRTAPs)
  "Q9BYR7", "Q9BYR6", "Q9BYR5", "Q9BYR4", "Q9BYR3", "Q9I7V9", "Q96D84",
  "Q07283", "P60331", "P60411", "P60413", "Q9BYQ5", "Q9BYP8", "Q8IUC1"
)

expr3 <- expr3[!(rownames(expr3) %in% keratin_ids), ]

# 5. Full interaction model
design <- model.matrix(~Cardiff_Batch_ID + AD_Status * condition, data = meta3)
colnames(design) <- make.names(colnames(design))

corfit <- duplicateCorrelation(expr3, design, block = meta3$donor)
fit <- lmFit(expr3, design, block = meta3$donor, correlation = corfit$consensus)
fit <- eBayes(fit,robust = TRUE)

# identify coefficients robustly
cn <- colnames(fit$coefficients)
int_coef <- grep("AD_Status.*[:\\.].*condition", cn, value = TRUE)
stopifnot(length(int_coef) == 1)

# extract tables
res_interaction <- topTable(fit, coef = int_coef, number = Inf, sort.by = "P")
res_AD_ctrl     <- topTable(fit, coef = "AD_StatusLOAD", number = Inf, sort.by = "P")
res_LPS_main    <- topTable(fit, coef = grep("^condition", cn, value = TRUE), number = Inf)

write.csv(res_interaction, file.path(outdir, "interaction_ADxLPS.csv"))
write.csv(res_AD_ctrl,     file.path(outdir, "LOAD_vs_CTRL_in_ctrl.csv"))
write.csv(res_LPS_main,    file.path(outdir, "LPS_main_effect.csv"))

# Contrasts 
cont <- makeContrasts(
  AD_in_ctrl = AD_StatusLOAD,
  AD_in_LPS  = AD_StatusLOAD + AD_StatusLOAD.conditionlps,
  Interaction = AD_StatusLOAD.conditionlps,
  # New Contrasts for the Plot
  LPS_Response_in_Healthy = conditionlps,
  LPS_Response_in_AD      = conditionlps + AD_StatusLOAD.conditionlps,
  
  levels = design
)

fit2 <- contrasts.fit(fit, cont)
fit2 <- eBayes(fit2,robust = TRUE)

res_AD_ctrl <- topTable(fit2, coef = "AD_in_ctrl", number = Inf)
res_AD_LPS  <- topTable(fit2, coef = "AD_in_LPS",  number = Inf)
res_int     <- topTable(fit2, coef = "Interaction", number = Inf)

write.csv(res_AD_LPS, file.path(outdir, "LOAD_vs_CTRL_in_LPS_contrast"))
write.csv(res_AD_ctrl,     file.path(outdir, "LOAD_vs_CTRL_in_ctrl_contrast.csv"))
write.csv(res_interaction,    file.path(outdir, "interaction_contrast.csv"))

sig_LPS <- res_AD_LPS %>%
  filter(adj.P.Val < 0.05)
sig_ctrl <- res_AD_ctrl %>%
  filter(adj.P.Val < 0.05)
sig_interaction <- res_int %>%
  filter(adj.P.Val < 0.05)

# PCA and volcano plots

# Calculate PCA on the averaged samples
pca_data <- expr3[order(apply(expr3, 1, var, na.rm=TRUE), decreasing=TRUE)[1:500], ]
pca <- prcomp(t(na.omit(pca_data)), scale. = TRUE)
pca_df <- as.data.frame(pca$x)
pca_df <- cbind(pca_df, meta3)
pca_df$Group <- interaction(pca_df$condition, pca_df$AD_Status)

# plot square aspect
png(paste0(outdir,"PCA.png"))
ggplot(pca_df, aes(x=PC1, y=PC2, color=Group, shape=Group)) +
  geom_point(size=3) +
  coord_fixed(ratio = 1) + 
  scale_color_manual(
    name = "Experimental Groups",
    values = c("ctrl.Control" = "#F8766D", "lps.Control" = "#00BFC4", 
               "ctrl.LOAD" = "#F8766D", "lps.LOAD" = "#00BFC4"),
    labels = c("Unstimulated LRWELL", "LPS Stimulated LRWELL", "Unstimulated HRLOAD", "LPS Stimulated HRLOAD")
  ) +
  scale_shape_manual(
    name = "Experimental Groups",
    values = c("ctrl.Control" = 16, "lps.Control" = 16, 
               "ctrl.LOAD" = 3, "lps.LOAD" = 3),
    labels = c("Unstimulated LRWELL", "LPS Stimulated LRWELL", "Unstimulated HRLOAD", "LPS Stimulated HRLOAD")
  ) +
  theme_minimal() +
  guides(color = guide_legend(override.aes = list(size = 4)))
dev.off()

plot_custom_volcano <- function(res_df, contrast_name, output_dir, show_labels = TRUE) {
  
  # Ensure numeric values
  res_df$P.Value   <- as.numeric(as.character(res_df$P.Value))
  res_df$adj.P.Val <- as.numeric(as.character(res_df$adj.P.Val))
  res_df$logFC     <- as.numeric(as.character(res_df$logFC))
  
  # Volcano significance colours
  keyvals <- rep('grey', nrow(res_df))
  is_up   <- which(!is.na(res_df$adj.P.Val) & res_df$adj.P.Val < 0.05 & res_df$logFC > 0)
  is_down <- which(!is.na(res_df$adj.P.Val) & res_df$adj.P.Val < 0.05 & res_df$logFC < 0)
  
  keyvals[is_up]   <- 'red3'
  keyvals[is_down] <- 'royalblue'
  key_names <- rep('NS', length(keyvals))
  key_names[is_up]   <- 'Signif. Up (adj.P.Val < 0.05)'
  key_names[is_down] <- 'Signif. Down (adj.P.Val < 0.05)'
  names(keyvals) <- key_names
  
  # Label display
  # Hide text without passing NULL to EnhancedVolcano
  current_lab_size <- if(show_labels) 4.5 else 0.0
  current_lab_col  <- if(show_labels) 'black' else 'transparent'
  selected_labels  <- if(show_labels) rownames(res_df)[keyvals != 'grey'] else NULL
  
  # Plot and save
  file_suffix <- if(show_labels) "_labeled.png" else "_unlabeled.png"
  file_path <- paste0(output_dir, "Volcano_", contrast_name, file_suffix)
  
  png(file_path, width = 900, height = 900, res = 130)
  
  p <- EnhancedVolcano(res_df,
                       lab = rownames(res_df), # Always provide to satisfy geom_text
                       x = 'logFC',
                       y = 'P.Value', 
                       selectLab = selected_labels,
                       title = contrast_name,
                       subtitle = if(show_labels) 'Labels: adj.P.Val < 0.05' else 'Labels: Hidden',
                       caption = '',
                       cutoffLineType = 'blank',
                       vline = NULL,
                       labSize = current_lab_size,
                       labCol = current_lab_col,
                       drawConnectors = show_labels,
                       widthConnectors = 0.4,
                       colCustom = keyvals,
                       colAlpha = 0.6,
                       pointSize = 1.5,
                       legendPosition = 'right') +
    theme(aspect.ratio = 1)
  
  print(p)
  dev.off()
  
  message(paste("Successfully saved:", file_path))
}

plot_custom_volcano(res_AD_LPS, "HRLOAD vs LRWELL LPS Stimulated", outdir, show_labels = TRUE)
plot_custom_volcano(res_LPS_main, "LPS Effect", outdir, show_labels = FALSE)

run_collapsed_fgsea <- function(res_table, outname, outdir) {
  # Map UniProt IDs to symbols
  res_table$Symbol <- AnnotationDbi::mapIds(org.Hs.eg.db,
                                            keys = sub(";.*", "", rownames(res_table)),
                                            column = "SYMBOL",
                                            keytype = "UNIPROT",
                                            multiVals = "first")
  
  # Rank genes
  ranks_df <- res_table %>%
    filter(!is.na(Symbol), !is.na(logFC), !is.na(P.Value)) %>%
    group_by(Symbol) %>%
    summarize(logFC = logFC[which.min(P.Value)],
              P.Value = min(P.Value)) %>%
    mutate(rank_score = sign(logFC) * -log10(P.Value + 1e-300))
  
  ranks <- setNames(ranks_df$rank_score, ranks_df$Symbol)
  ranks <- sort(ranks, decreasing = TRUE)
  
  # Reactome pathways
  m_df <- msigdbr::msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
  pathways <- split(m_df$gene_symbol, m_df$gs_name)
  
  # Enrichment
  set.seed(42) 
  fg_full <- fgsea::fgsea(pathways = pathways, 
                          stats = ranks, 
                          minSize = 10,   
                          maxSize = 500)
  
  # Collapse redundant pathways at nominal p < 0.05
  sig_fg <- fg_full[fg_full$pval < 0.05, ]
  
  if(nrow(sig_fg) > 0) {
    set.seed(42) # Added seed for reproducible collapse
    collapsed_pathways <- fgsea::collapsePathways(fgseaRes = sig_fg, 
                                                  pathways = pathways, 
                                                  stats = ranks)
    fg_collapsed <- fg_full[pathway %in% collapsed_pathways$mainPathways]
  } else {
    message("No significant pathways found to collapse.")
    fg_collapsed <- fg_full
  }
  
  # Save nominally significant pathways
  fg_final <- fg_collapsed[fg_collapsed$pval < 0.05, ]
  fg_final$leadingEdge <- sapply(fg_final$leadingEdge, function(x) paste(x, collapse = ";"))
  fg_final <- fg_final[order(fg_final$pval), ]
  
  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  write.csv(fg_final, file = file.path(outdir, outname), row.names = FALSE)
  
  return(fg_final)
}

fg_unstim_collapsed   <- run_collapsed_fgsea(res_AD_ctrl, "FGSEA_unstim_Collapsed.csv", outdir)
fg_LPS_collapsed   <- run_collapsed_fgsea(res_AD_LPS, "FGSEA_LPS_Collapsed.csv", outdir)

plot_gsea_results <- function(fgsea_res, title_str, outname, outdir, num_pathways = 10) {
  
  # Top pathways ranked by nominal p-value
  plot_data <- fgsea_res %>%
    arrange(pval) %>%
    slice_head(n = num_pathways) %>%
    mutate(
      pathway = gsub("REACTOME_", "", pathway),
      pathway = gsub("_", " ", pathway),
      sig_star = if_else(padj < 0.05, "*", ""),
      Direction = if_else(NES > 0, "Up", "Down")
    )
  
  # Plot
  p <- ggplot(plot_data, aes(x = reorder(pathway, NES), y = NES, fill = Direction)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sig_star, 
                  hjust = if_else(NES > 0, -0.2, 1.2)), 
              vjust = 0.7, 
              size = 8, 
              fontface = "bold") +
    coord_flip() +
    scale_fill_manual(values = c("Up" = "#E41A1C", "Down" = "#377EB8")) +
      scale_y_continuous(expand = expansion(mult = c(0.2, 0.2))) +
    labs(
      title = title_str,
      subtitle = paste("Top", num_pathways, "Pathways by p-value | * padj < 0.05"),
      x = NULL,
      y = "Normalized Enrichment Score (NES)"
    ) +
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 9, face = "bold"),
      panel.grid.minor = element_blank(),
      legend.position = "bottom"
    )
  
  # Save
  if(!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
  file_path <- file.path(outdir, outname)
  
  ggsave(file_path, plot = p, width = 8, height = 6, dpi = 150)
  
  message(paste("GSEA plot saved to:", file_path))
  return(p)
}

plot_gsea_results(fg_LPS_collapsed, "HRLOAD vs LRWELL (LPS Stimulated)", "GSEA_LPS_Plot.png", outdir)

plot_gsea_results(fg_unstim_collapsed, "HRLOAD vs LRWELL (Unstimulated)", "GSEA_Unstimulated_Plot.png", outdir)


# Diverse Reactome enrichment plots

plot_diverse_gsea <- function(df, title_str, outname, outdir, overlap_cutoff = 0.4) {
  
  if (nrow(df) == 0) return(message("Empty dataframe for ", title_str))
  
  # Select pathways with limited leading-edge gene overlap
  df <- df[order(df$pval), ]
  keep_indices <- c(1)
  get_genes <- function(idx) unlist(strsplit(as.character(df$leadingEdge[idx]), ";"))
  
  if (nrow(df) > 1) {
    for (i in 2:nrow(df)) {
      current_genes <- get_genes(i)
      is_redundant <- FALSE
      for (j in keep_indices) {
        jaccard <- length(intersect(current_genes, get_genes(j))) / 
          length(union(current_genes, get_genes(j)))
        if (jaccard > overlap_cutoff) { is_redundant <- TRUE; break }
      }
      if (!is_redundant) keep_indices <- c(keep_indices, i)
      if (length(keep_indices) == 10) break
    }
  }
  
  plot_df <- df[keep_indices, ] %>%
    mutate(
      short_name = gsub("REACTOME_", "", pathway) %>% gsub("_", " ", .),
      is_sig = (padj < 0.05),
      sig_label = if_else(is_sig, "*", ""),
      Direction = if_else(NES > 0, "Up", "Down")
    )
  
  # Plot
  p <- ggplot(plot_df, aes(x = NES, y = reorder(short_name, NES), fill = Direction)) +
    geom_col(alpha = 0.8) +
    geom_text(aes(label = sig_label, 
                  hjust = if_else(NES > 0, -0.3, 1.3)), 
              vjust = 0.8, size = 8, fontface = "bold") +
    scale_fill_manual(values = c("Up" = "#e41a1c", "Down" = "#377eb8")) +
    scale_x_continuous(expand = expansion(mult = c(0.3, 0.4))) +
    theme_minimal() +
    labs(title = title_str,
         subtitle = paste0("Max ", overlap_cutoff*100, "% gene overlap | * = padj < 0.05"),
         x = "Normalized Enrichment Score (NES)", y = NULL) +
    theme(
      axis.text.y = element_text(size = 9, face = "bold"),
      legend.position = "bottom",
      plot.margin = margin(t = 10, r = 100, b = 10, l = 10)
    )
  
  ggsave(file.path(outdir, outname), p, width = 11, height = 7, dpi = 150)
}

for(s in c("LPS","unstim")) {
  fname <- file.path(outdir, paste0("FGSEA_", s, "_Collapsed.csv"))
  
  if(file.exists(fname)) {
    raw_res <- read.csv(fname, stringsAsFactors = FALSE)
    
    plot_diverse_gsea(
      df = raw_res, 
      title_str = paste("Diverse Top 10:", toupper(s)),
      outname = paste0("BarChart_Diverse_Final_", s, ".png"),
      outdir = outdir,
      overlap_cutoff = 0.4
    )
  }
}

save.image("full_proteomics_with_sample_exclusions.Rdata")
