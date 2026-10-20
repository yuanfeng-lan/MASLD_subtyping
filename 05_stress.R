rm(list = ls()); gc()

## -------- Paths & directories --------
ORIGINAL_DIR <- ""                        # working directory (set as needed)
output       <- file.path(ORIGINAL_DIR, "05_stress")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
setwd(ORIGINAL_DIR)

# Purpose: GSVA scoring of cellular-stress gene sets (oxidative stress, hypoxia,
# chronic inflammation, ER stress, intracellular/extracellular pH elevation),
# then characterize stress states across subtypes and their relation to metabolism.
# k = 4.

## -------- Libraries --------
library(GSVA)
library(GSEABase)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(grid)

## -------- Settings --------
labels   <- "05_stress_"
cols     <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")
gse_list <- c("GSE193066", "GSE213621", "GSE83452", "GSE164760")

# Metabolic / ECM gene sets used for the stress-vs-metabolism correlation heatmaps
metab_sets <- c("GOBP_AMINO_ACID_METABOLIC_PROCESS", "GOBP_AMINO_ACID_CATABOLIC_PROCESS",
                "GOBP_FATTY_ACID_METABOLIC_PROCESS", "GOBP_FATTY_ACID_CATABOLIC_PROCESS",
                "GOBP_GOLGI_VESICLE_TRANSPORT",
                "GOBP_POSITIVE_REGULATION_OF_MRNA_CATABOLIC_PROCESS",
                "GOBP_DNA_CATABOLIC_PROCESS", "GOBP_REGULATION_OF_DNA_METABOLIC_PROCESS")

# Stress axes and metabolic axes used in the scatter plots
x_list <- c("HYPOXIA", "CHRONIC_INFLAMMATORY", "OXIDATIVE_STRESS", "ER_STRESS")
y_list <- c("GOBP_FATTY_ACID_CATABOLIC_PROCESS", "GOBP_AMINO_ACID_CATABOLIC_PROCESS",
            "GOBP_FATTY_ACID_METABOLIC_PROCESS", "GOBP_AMINO_ACID_METABOLIC_PROCESS",
            "GOBP_GOLGI_VESICLE_TRANSPORT", "GOBP_REGULATION_OF_GLYCOGEN_METABOLIC_PROCESS",
            "GOBP_REGULATION_OF_GLYCOGEN_BIOSYNTHETIC_PROCESS",
            "GOBP_POSITIVE_REGULATION_OF_GLYCOGEN_BIOSYNTHETIC_PROCESS",
            "GOBP_POSITIVE_REGULATION_OF_GLYCOPROTEIN_METABOLIC_PROCESS")

## ============================================================================
## Helper: read a single-line GMT and collapse to a unique gene vector
## ============================================================================
read_gmt_genes <- function(path) {
  gmt <- read.csv(path, header = FALSE, sep = "\t")
  gmt <- gmt[-1, ]
  gene_set <- gmt[, -2]
  rownames(gene_set) <- NULL
  gene_set <- gene_set %>% column_to_rownames("V1") %>% t()
  genes <- unique(as.character(unlist(gene_set)))
  genes[nchar(genes) != 0]
}

## ============================================================================
## Helper: GSVA over a gene-set list, return a data.frame
## ============================================================================
run_gsva <- function(expr, gene_set) {
  param <- gsvaParam(as.matrix(expr), gene_set)
  as.data.frame(gsva(param, verbose = TRUE))
}

## ============================================================================
## Helper: per-subtype Spearman correlation between stress and metabolic sets
## ============================================================================
cor_heatmap <- function(stress_data, meta_data, group, out_dir, prefix = "") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  for (CS in unique(group$clust)) {
    sample  <- rownames(group[group$clust == CS, ])
    data1   <- stress_data[sample, ]
    data2   <- meta_data[sample, ]
    data2   <- data2[rownames(data1), ]

    p01 <- colnames(data1); p02 <- colnames(data2)
    r_mat <- matrix(NA_real_, nrow = length(p01), ncol = length(p02),
                    dimnames = list(p01, p02))
    p_mat <- r_mat

    for (i in seq_along(p01)) for (j in seq_along(p02)) {
      x <- data1[[p01[i]]]; y <- data2[[p02[j]]]
      ok <- complete.cases(x, y)
      if (sum(ok) < 3 || sd(x[ok]) == 0 || sd(y[ok]) == 0) next
      ct <- cor.test(x[ok], y[ok], method = "spearman")
      r_mat[i, j] <- as.numeric(ct$estimate)
      p_mat[i, j] <- ct$p.value
    }

    write.csv(as.data.frame(r_mat), file.path(out_dir, paste0(prefix, CS, "_pathway_cor_r.csv")))
    write.csv(as.data.frame(p_mat), file.path(out_dir, paste0(prefix, CS, "_pathway_cor_p.csv")))

    r_long <- as.data.frame(as.table(r_mat)); colnames(r_long) <- c("pathway1", "pathway2", "r")
    p_long <- as.data.frame(as.table(p_mat)); colnames(p_long) <- c("pathway1", "pathway2", "p")
    plot_df <- left_join(r_long, p_long, by = c("pathway1", "pathway2"))
    plot_df$star <- with(plot_df, ifelse(is.na(p),
                                         "",
                                         ifelse(p <= 0.001, "***",
                                                ifelse(p <= 0.01, "**",
                                                       ifelse(p <= 0.05, "*", "")))))
    plot_df$pathway1 <- factor(plot_df$pathway1, levels = p01)
    plot_df$pathway2 <- factor(plot_df$pathway2, levels = p02)

    p_heat <- ggplot(plot_df, aes(x = pathway1, y = pathway2, fill = r)) +
      geom_tile(color = "white", size = 0.2) +
      scale_fill_gradient2(low = "#3C5488", mid = "white", high = "#E64B35",
                           midpoint = 0, limits = c(-1, 1), na.value = "grey90") +
      geom_text(aes(label = star), color = "black", size = 3) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
            axis.text.y = element_text(hjust = 1, vjust = 1, face = "bold"),
            axis.title = element_blank(),
            panel.grid = element_blank()) +
      labs(fill = "Spearman r", title = CS)
    ggsave(file.path(out_dir, paste0(prefix, CS, "pathway_cor_heatmap.png")),
           p_heat, width = 8, height = 6, dpi = 300)
    ggsave(file.path(out_dir, paste0(prefix, CS, "pathway_cor_heatmap.pdf")),
           p_heat, width = 8, height = 6)
  }
}

## ============================================================================
## Helper: stress-score boxplot by subtype (Kruskal-Wallis p)
## ============================================================================
plot_stress_box <- function(df, pathway, group_col = "group",
                            cols, out_dir, prefix = "", subtitle = NULL,
                            width = 4, height = 4) {
  fml  <- as.formula(paste(pathway, "~", group_col))
  pval <- compare_means(fml, data = df, method = "kruskal.test")$p.format
  ymax <- max(df[[pathway]], na.rm = TRUE)

  p <- ggplot(df, aes(x = .data[[group_col]], y = .data[[pathway]],
                      fill = .data[[group_col]], colour = .data[[group_col]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5,
                 colour = "black", size = 0.4) +
    geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    theme_minimal(base_size = 15) +
    labs(x = "", y = "Pathway GSVA Score", title = pathway, subtitle = subtitle) +
    theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.title.x    = element_text(size = 16, face = "bold"),
          axis.title.y    = element_text(size = 14),
          axis.text.x     = element_text(size = 14, face = "bold"),
          axis.text.y     = element_text(size = 14),
          legend.position = "none",
          plot.title      = element_text(hjust = 0.5, size = 16, face = "bold")) +
    coord_cartesian(clip = "off") +
    annotation_custom(
      grob = textGrob(label = paste0("p = ", pval), gp = gpar(fontsize = 15), hjust = 0.5),
      xmin = 1.1, xmax = 1.5,
      ymin = ymax * 0.85, ymax = ymax * 1.05)

  ggsave(file.path(out_dir, paste0(prefix, pathway, "_plot.png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(prefix, pathway, "_plot.pdf")), p, width = width, height = height)
}

## ============================================================================
## Helper: stress-vs-metabolism scatter with Pearson correlation
## ============================================================================
plot_stress_metab_scatter <- function(stress_data, meta_data, x_list, y_list,
                                      out_dir, prefix = "", subtitle = NULL,
                                      width = 4, height = 4) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  for (i in x_list) for (m in y_list) {
    x <- stress_data[[i]]; y <- meta_data[[m]]
    data <- na.omit(data.frame(x = x, y = y))

    p <- ggplot(data, aes(x = x, y = y)) +
      geom_point(color = "#4DBBD5", size = 2, alpha = 0.9) +
      geom_smooth(method = "lm", color = "#E64B35", se = TRUE, linewidth = 0.8) +
      stat_cor(method = "pearson",
               label.x = min(data$x), label.y = max(data$y) * 0.9, size = 5) +
      labs(x = gsub("_", " ", i), y = gsub("_", " ", m), subtitle = subtitle) +
      theme_minimal(base_size = 15) +
      theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
            axis.title.x    = element_text(size = 12, face = "bold"),
            axis.title.y    = element_text(size = 12, face = "bold"),
            axis.text.x     = element_text(size = 14, face = "bold"),
            axis.text.y     = element_text(size = 14),
            legend.position = "none",
            plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle   = element_text(hjust = 0.5, size = 14)) +
      coord_cartesian(clip = "off")
    ggsave(file.path(out_dir, paste0(prefix, i, "_vs_", m, "_plot.png")), p, width = width, height = height, dpi = 300)
    ggsave(file.path(out_dir, paste0(prefix, i, "_vs_", m, "_plot.pdf")), p, width = width, height = height)
  }
}

## ############################################################################
## 1. Build stress gene sets from GMT files
## ############################################################################
gmt_root <- file.path(output, "gmt")
gene_set <- list()

gene_set[["OXIDATIVE_STRESS"]] <- read_gmt_genes(
  file.path(gmt_root, "oxidative_stress_geneset", "genesets.v2026.1.Hs.gmt"))
gene_set[["HYPOXIA"]] <- read_gmt_genes(
  file.path(gmt_root, "hypoxia_geneset", "genesets.v2026.1.Hs.gmt"))
gene_set[["CHRONIC_INFLAMMATORY"]] <- read_gmt_genes(
  file.path(gmt_root, "chronic_inflammation_geneset", "genesets.v2026.1.Hs.gmt"))
gene_set[["ER_STRESS"]] <- read_gmt_genes(
  file.path(gmt_root, "er_stress_geneset", "genesets.v2026.1.Hs.gmt"))
gene_set[["INTRACELLULAR_PH_ELEVATION"]] <- read_gmt_genes(
  file.path(gmt_root, "stress_geneset", "5-pHi.gmt"))
gene_set[["EXTRACELLULAR_PH_ELEVATION"]] <- read_gmt_genes(
  file.path(gmt_root, "stress_geneset", "4-pHe.gmt"))

save(gene_set, file = file.path(output, paste0(labels, "Stess_gmt.Rdata")))

max_len <- max(sapply(gene_set, length))
combined_df <- data.frame(
  OXIDATIVE_STRESS            = c(gene_set$OXIDATIVE_STRESS,            rep(NA, max_len - length(gene_set$OXIDATIVE_STRESS))),
  HYPOXIA                     = c(gene_set$HYPOXIA,                     rep(NA, max_len - length(gene_set$HYPOXIA))),
  CHRONIC_INFLAMMATORY        = c(gene_set$CHRONIC_INFLAMMATORY,        rep(NA, max_len - length(gene_set$CHRONIC_INFLAMMATORY))),
  ER_STRESS                   = c(gene_set$ER_STRESS,                   rep(NA, max_len - length(gene_set$ER_STRESS))),
  INTRACELLULAR_PH_ELEVATION  = c(gene_set$INTRACELLULAR_PH_ELEVATION,  rep(NA, max_len - length(gene_set$INTRACELLULAR_PH_ELEVATION))),
  EXTRACELLULAR_PH_ELEVATION  = c(gene_set$EXTRACELLULAR_PH_ELEVATION,  rep(NA, max_len - length(gene_set$EXTRACELLULAR_PH_ELEVATION))))
write.csv(combined_df, file = file.path(output, paste0(labels, "stress_gene_sets.csv")))

## ############################################################################
## 2. GSVA on the discovery cohort
## ############################################################################
load(file.path(ORIGINAL_DIR, "00_rawdata_2", "exp_merge_combat_nonneg.rData"))
exp_MASLD <- as.data.frame(res)

group <- read.csv(file.path(ORIGINAL_DIR, "01_MOVICS", "01_MOVICS_cluster_results.csv"), row.names = 1)
group$clust <- paste0("CS", group$clust)

gsva_results <- run_gsva(exp_MASLD, gene_set)
save(gsva_results, file = file.path(output, paste0(labels, "stress_gsva_results.Rdata")))
write.csv(gsva_results, file.path(output, paste0(labels, "stress_gsva_results.csv")))

## ############################################################################
## 3. WGCNA module vs. stress-trait correlation heatmap
## ############################################################################
trait_data <- as.data.frame(t(gsva_results))
load(file.path(ORIGINAL_DIR, "03_WGCNA_GSVA", "03_WGCNA_GSVA_WGCNA_results.RData"))
load(file.path(ORIGINAL_DIR, "02_GSVA_k4",  "02_GSVA_k4_metabolic_gave.Rdata"))

datExpr <- as.data.frame(t(metabolic_gave))
datExpr <- as.data.frame(lapply(datExpr, as.numeric))
rownames(datExpr) <- colnames(metabolic_gave)

moduleTraitCor    <- cor(MEs, trait_data, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

draw_module_trait <- function() {
  par(mar = c(5, 10, 4, 2))
  labeledHeatmap(Matrix      = moduleTraitCor,
                 xLabels     = names(trait_data),
                 yLabels     = names(MEs),
                 ySymbols    = names(MEs),
                 colorLabels = FALSE,
                 colors      = blueWhiteRed(50),
                 textMatrix  = textMatrix,
                 setStdMargins = FALSE,
                 cex.text    = 0.5,
                 zlim        = c(-1, 1),
                 main        = "Module-Cluster relationships")
}
png(file.path(output, paste0(labels, "WGCNA_Module-trait relationships_plot.png")),
    width = 5, height = 8, res = 300, units = "in"); draw_module_trait(); dev.off()
pdf(file.path(output, paste0(labels, "WGCNA_Module-trait relationships_plot.pdf")),
    width = 5, height = 8); draw_module_trait(); dev.off()

## ############################################################################
## 4. Stress-metabolism correlation heatmaps (discovery)
## ############################################################################
load(file.path(ORIGINAL_DIR, "02_GSVA_k4", "02_GSVA_k4_gsva_results.Rdata"))
meta_data <- as.data.frame(t(gsva_results))[, metab_sets]

cor_heatmap(as.data.frame(t(gsva_results)), meta_data, group,
            out_dir = file.path(output, "Cor_test"))

## ############################################################################
## 5. GSVA on external validation cohorts
## ############################################################################
for (GSE in gse_list) {
  exp   <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(GSE, "_exp.csv")), row.names = 1)
  out_dir <- file.path(output, GSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  gsva_results <- run_gsva(exp, gene_set)
  save(gsva_results, file = file.path(out_dir, paste0(labels, GSE, "_stress_gsva_results.Rdata")))
  write.csv(gsva_results, file.path(out_dir, paste0(labels, GSE, "_stress_gsva_results.csv")))
}

## ############################################################################
## 6. Stress-metabolism correlation heatmaps (external cohorts)
## ############################################################################
for (GSE in gse_list) {
  out_dir <- file.path(output, "Cor_test", GSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  load(file.path(output, GSE, paste0(labels, GSE, "_stress_gsva_results.Rdata")))
  stress_data <- as.data.frame(t(gsva_results))

  group  <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(GSE, "_group.csv")), row.names = 1)
  cluster_group <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(labels, GSE, "cluster_res.csv")), row.names = 1)
  cluster_group$clust <- paste0("CS", cluster_group$clust)
  cluster_group <- cluster_group[rownames(group), ]
  group$clust   <- cluster_group$clust

  load(file.path(ORIGINAL_DIR, "02_GSVA_k4", GSE, paste0("02_GSVA_k4_", GSE, "_gobp_gsva_results.Rdata")))
  meta_data <- as.data.frame(t(gsva_results))[, metab_sets]

  cor_heatmap(stress_data, meta_data, group, out_dir, prefix = GSE)
}

## ############################################################################
## 7. Stress-level boxplots by subtype (discovery)
## ############################################################################
load(file.path(output, paste0(labels, "stress_gsva_results.Rdata")))
stress_data <- as.data.frame(t(gsva_results))

group <- read.csv(file.path(ORIGINAL_DIR, "01_MOVICS", "01_MOVICS_cluster_results.csv"), row.names = 1)
group$clust <- paste0("CS", group$clust)
group <- group[rownames(stress_data), ]

data_df <- stress_data
data_df$group <- group$clust

out_dir <- file.path(output, "stress_test")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
for (pathway in colnames(stress_data)) {
  plot_stress_box(data_df, pathway, cols = cols, out_dir = out_dir, prefix = labels)
}

## ############################################################################
## 8. Stress-level boxplots by subtype (external cohorts)
## ############################################################################
for (GSE in gse_list) {
  out_dir <- file.path(output, "stress_test", GSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  load(file.path(output, GSE, paste0(labels, GSE, "_stress_gsva_results.Rdata")))
  stress_data <- as.data.frame(t(gsva_results))

  group <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(GSE, "_group.csv")), row.names = 1)
  cluster_group <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(labels, GSE, "cluster_res.csv")), row.names = 1)
  cluster_group$clust <- paste0("CS", cluster_group$clust)
  cluster_group <- cluster_group[rownames(group), ]
  group$clust   <- cluster_group$clust

  data_df <- stress_data
  data_df$group <- group$clust
  for (pathway in colnames(stress_data)) {
    plot_stress_box(data_df, pathway, cols = cols, out_dir = out_dir,
                    prefix = GSE, subtitle = GSE)
  }
}

## ############################################################################
## 9. Stress-vs-metabolism scatter plots (discovery)
## ############################################################################
load(file.path(output, paste0(labels, "stress_gsva_results.Rdata")))
stress_data <- as.data.frame(t(gsva_results))
load(file.path(ORIGINAL_DIR, "02_GSVA_k4", "02_GSVA_k4_gsva_results.Rdata"))
meta_data <- as.data.frame(t(gsva_results))[, metab_sets]

plot_stress_metab_scatter(stress_data, meta_data, x_list, y_list,
                          out_dir = file.path(output, "Cor_test_2"))

## ############################################################################
## 10. Stress-vs-metabolism scatter plots (external cohorts)
## ############################################################################
ext_metab_sets <- c(metab_sets,
                    "GOBP_EXTRACELLULAR_MATRIX_CONSTITUENT_SECRETION",
                    "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY",
                    "GOBP_POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_DISASSEMBLY",
                    "GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_CONSTITUENT_SECRETION",
                    "GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION",
                    "GOBP_REGULATION_OF_GLYCOGEN_METABOLIC_PROCESS",
                    "GOBP_REGULATION_OF_GLYCOGEN_BIOSYNTHETIC_PROCESS",
                    "GOBP_POSITIVE_REGULATION_OF_GLYCOGEN_BIOSYNTHETIC_PROCESS",
                    "GOBP_POSITIVE_REGULATION_OF_GLYCOPROTEIN_METABOLIC_PROCESS")

for (GSE in gse_list) {
  out_dir <- file.path(output, "Cor_test_2", GSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  load(file.path(output, GSE, paste0(labels, GSE, "_stress_gsva_results.Rdata")))
  stress_data <- as.data.frame(t(gsva_results))

  load(file.path(ORIGINAL_DIR, "02_GSVA_k4", GSE, paste0("02_GSVA_k4_", GSE, "_gobp_gsva_results.Rdata")))
  meta_data <- as.data.frame(t(gsva_results))[, ext_metab_sets]

  plot_stress_metab_scatter(stress_data, meta_data, x_list, y_list,
                            out_dir = out_dir, prefix = GSE, subtitle = GSE)
}

## ############################################################################
## 11. Radar plots of cellular stress by subtype
## ############################################################################
load(file.path(output, paste0(labels, "stress_gsva_results.Rdata")))
group <- read.csv(file.path(ORIGINAL_DIR, "01_MOVICS", "01_MOVICS_cluster_results.csv"), row.names = 1)
group$clust <- paste0("CS", group$clust)

group_levels  <- c("CS1", "CS2", "CS3", "CS4")
stress_order  <- rownames(gsva_results)
output_radar  <- file.path(output, "radar")
if (!dir.exists(output_radar)) dir.create(output_radar, recursive = TRUE)

group_df <- data.frame(sample = rownames(group), group = group$clust, stringsAsFactors = FALSE)
group_df <- group_df %>%
  dplyr::distinct() %>%
  dplyr::mutate(group = factor(group, levels = group_levels)) %>%
  dplyr::arrange(group, sample)
rownames(group_df) <- group_df$sample

common_samples <- intersect(group_df$sample, colnames(gsva_results))
if (length(common_samples) == 0) stop("No matched samples between group file and gsva_results.")
if (length(common_samples) < nrow(group_df)) warning("Some group samples are missing in gsva_results and will be dropped.")
group_df <- group_df[common_samples, , drop = FALSE]

gsva_long <- as.data.frame(gsva_results) %>%
  tibble::rownames_to_column("Function") %>%
  tidyr::pivot_longer(cols = -Function, names_to = "sample", values_to = "GSVA") %>%
  dplyr::filter(sample %in% common_samples) %>%
  dplyr::left_join(group_df %>% dplyr::select(sample, group), by = "sample")
if (any(is.na(gsva_long$group))) stop("Some samples do not match the group file. Please check sample names.")

gsva_cat_sample <- gsva_long %>%
  dplyr::group_by(sample, group, Function) %>%
  dplyr::summarise(Function_GSVA = mean(GSVA, na.rm = TRUE), .groups = "drop")

gsva_cat_group <- gsva_cat_sample %>%
  dplyr::group_by(group, Function) %>%
  dplyr::summarise(mean_GSVA = mean(Function_GSVA, na.rm = TRUE), .groups = "drop")

baseline_df <- gsva_cat_sample %>%
  dplyr::group_by(Function) %>%
  dplyr::summarise(baseline_GSVA = mean(Function_GSVA, na.rm = TRUE), .groups = "drop")

pairwise_res <- bind_rows(lapply(stress_order, function(fn) {
  df_f <- gsva_cat_sample %>%
    dplyr::filter(Function == fn) %>%
    dplyr::mutate(group = factor(group, levels = group_levels))
  bind_rows(lapply(group_levels, function(g) {
    x <- df_f %>% dplyr::filter(group == g) %>% dplyr::pull(Function_GSVA)
    y <- df_f %>% dplyr::filter(group != g) %>% dplyr::pull(Function_GSVA)
    if (length(x) < 2 || length(y) < 2) return(NULL)
    p <- tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value, error = function(e) NA_real_)
    tibble::tibble(Function = fn, group = g, p_value = p,
                   median_group = median(x, na.rm = TRUE),
                   median_other = median(y, na.rm = TRUE),
                   direction    = ifelse(median(x, na.rm = TRUE) > median(y, na.rm = TRUE), "Higher", "Lower"))
  }))
})) %>%
  dplyr::mutate(fdr  = p.adjust(p_value, method = "BH"),
                star = dplyr::case_when(fdr < 0.001 ~ "***",
                                        fdr < 0.01  ~ "**",
                                        fdr < 0.05  ~ "*",
                                        TRUE ~ ""))

write.csv(pairwise_res, file.path(output_radar, "Cell_stress_one_vs_rest_significance.csv"), row.names = FALSE)

y_rng <- range(c(gsva_cat_group$mean_GSVA, baseline_df$baseline_GSVA), na.rm = TRUE)
y_pad <- diff(y_rng) * 0.10; if (is.na(y_pad) || y_pad == 0) y_pad <- 0.1
y_lim <- c(y_rng[1] - y_pad, y_rng[2] + y_pad)

make_radar_plot <- function(group_name, plot_df, pairwise_res, baseline_df, y_lim, output_dir) {
  df_g <- plot_df %>%
    dplyr::filter(group == group_name) %>%
    dplyr::select(Function, mean_GSVA) %>%
    tidyr::complete(Function = factor(stress_order, levels = stress_order),
                    fill = list(mean_GSVA = 0)) %>%
    dplyr::mutate(Function = factor(Function, levels = stress_order),
                  x = seq_along(stress_order)) %>%
    dplyr::arrange(Function)
  df_closed <- dplyr::bind_rows(df_g, dplyr::mutate(df_g[1, ], x = nrow(df_g) + 1))

  baseline_plot_df <- baseline_df %>%
    dplyr::select(Function, baseline_GSVA) %>%
    tidyr::complete(Function = factor(stress_order, levels = stress_order),
                    fill = list(baseline_GSVA = 0)) %>%
    dplyr::mutate(Function = factor(Function, levels = stress_order),
                  x = seq_along(stress_order)) %>%
    dplyr::arrange(Function)
  baseline_closed <- dplyr::bind_rows(baseline_plot_df,
                                      dplyr::mutate(baseline_plot_df[1, ], x = nrow(baseline_plot_df) + 1))

  sig_df <- pairwise_res %>%
    dplyr::filter(group == group_name, star != "") %>%
    dplyr::distinct(Function, .keep_all = TRUE) %>%
    dplyr::mutate(x = match(Function, stress_order),
                  y = y_lim[2] + 0.06 * diff(y_lim))

  p <- ggplot() +
    geom_path(data = baseline_closed, aes(x = x, y = baseline_GSVA, group = 1),
              color = "grey40", linetype = "dashed", linewidth = 1.0) +
    geom_path(data = df_closed, aes(x = x, y = mean_GSVA, group = 1),
              color = "#2C7FB8", linewidth = 1.0) +
    geom_point(data = df_g, aes(x = x, y = mean_GSVA), color = "#2C7FB8", size = 2.2) +
    geom_text(data = sig_df, aes(x = x, y = y, label = star),
              inherit.aes = FALSE, color = "red", size = 6, fontface = "bold") +
    coord_polar(start = -pi / 2, clip = "off") +
    scale_x_continuous(breaks = seq_along(stress_order), labels = stress_order,
                       limits = c(1, nrow(df_g) + 1)) +
    scale_y_continuous(limits = c(y_lim[1], y_lim[2] + 0.12 * diff(y_lim)),
                       breaks = pretty(y_lim, n = 5)) +
    labs(title = paste0(group_name, " Cell Stress GSVA Radar Plot"),
         x = NULL, y = "Raw GSVA") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_line(color = "grey85"),
          axis.text.x  = element_text(size = 9, color = "black"),
          axis.text.y  = element_text(size = 9, color = "grey30"),
          axis.title.y = element_text(size = 11),
          plot.title   = element_text(hjust = 0.5, face = "bold"),
          plot.margin  = margin(15, 20, 15, 20))

  ggsave(file.path(output_dir, paste0(group_name, "_cell_stress_radar.png")), p, width = 10, height = 8, dpi = 300)
  ggsave(file.path(output_dir, paste0(group_name, "_cell_stress_radar.pdf")), p, width = 10, height = 8)
  p
}

plots_list <- list()
for (g in group_levels) {
  plots_list[[g]] <- make_radar_plot(g, gsva_cat_group, pairwise_res, baseline_df, y_lim, output_radar)
}
write.csv(gsva_cat_sample, file.path(output_radar, "cell_stress_each_sample_raw.csv"), row.names = FALSE)
write.csv(gsva_cat_group,  file.path(output_radar, "cell_stress_group_mean_raw.csv"), row.names = FALSE)
write.csv(baseline_df,     file.path(output_radar, "cell_stress_baseline_raw.csv"), row.names = FALSE)
