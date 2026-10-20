rm(list = ls()); gc()

## -------- Paths & directories --------
ORIGINAL_DIR <- ""                        # working directory (set as needed)
output       <- file.path(ORIGINAL_DIR, "06_CIBERSORT")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
setwd(ORIGINAL_DIR)

# Purpose: immune-infiltration analysis across subtypes (CIBERSORT, ESTIMATE,
# MCPcounter), with particular focus on whether CS3 shows higher immune activation.
# k = 4.

## -------- Libraries --------
library(tidyverse)
library(e1071)
library(parallel)
library(preprocessCore)
library(reshape2)
library(aplot)
library(psych)
library(corrplot)
library(ggpubr)
library(grid)
library(rstatix)
library(IOBR)
library(dplyr)
library(tibble)

## -------- Settings --------
labels       <- "06_CIBERSORT_"
gse_list     <- c("GSE193066", "GSE213621", "GSE83452", "GSE164760")
disc_cols    <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")

# CIBERSORT implementation and signature matrix
source(file.path(output, "CIBERSORT.R"))
sig_matrix <- file.path(output, "LM22.txt")

## ============================================================================
## Helper: run CIBERSORT on an expression matrix, return the 22-column result
## ============================================================================
run_cibersort <- function(exp, out_dir, tag) {
  write.table(exp, file.path(out_dir, paste0(tag, "CIBERSORT_exp.csv")),
              sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
  mixture_file <- file.path(out_dir, paste0(tag, "CIBERSORT_exp.csv"))
  res <- CIBERSORT(sig_matrix, mixture_file, perm = 100, QN = TRUE)
  res <- res[, 1:22]
  ciber.res <- res[, colSums(res) > 0]
  write.csv(ciber.res, file.path(out_dir, paste0(labels, tag, "CIBERSORT_ciber.res.csv")))
  ciber.res
}

## ============================================================================
## Helper: stacked-bar composition plot + subtype annotation strip
## ============================================================================
plot_composition <- function(ciber.res, group, out_dir, labels, tag, title = NULL) {
  group <- group[order(group$clust), ]
  group$x <- 1
  group$y <- group$clust
  group$sample <- rownames(group)

  p2 <- ggplot(group, aes(x = x, y = sample)) +
    geom_tile(aes(fill = clust)) +
    scale_x_continuous(expand = c(0, 0)) +
    theme(panel.background = element_blank(),
          axis.ticks       = element_blank(),
          axis.title       = element_blank(),
          axis.text.x      = element_blank(),
          axis.text.y      = element_blank(),
          legend.position  = "left",
          legend.title     = element_blank()) +
    scale_fill_manual(values = c("#CC88B0", "#998DB7", "#DBE0ED", "#87B5B2", "#F4CEB4")) +
    coord_flip()

  ciber.res <- ciber.res[rownames(group), ]
  mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7)

  ciber_res_long <- melt(ciber.res)
  colnames(ciber_res_long) <- c("sample", "cell type", "Relative percentage")
  ciber_res_long$sample <- factor(ciber_res_long$sample, levels = group$sample)

  p1 <- ggplot(ciber_res_long, aes(x = sample, y = `Relative percentage`)) +
    geom_bar(stat = "identity", position = "stack", aes(fill = `cell type`)) +
    labs(x = NULL, title = title) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = mycol, name = NULL) +
    theme_bw() +
    theme(axis.text.x  = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "right")

  p <- insert_bottom(p1, p2, height = 0.05)
  ggsave(file.path(out_dir, paste0(labels, tag, "result01.pdf")), p, width = 8, height = 5)
  ggsave(file.path(out_dir, paste0(labels, tag, "result01.png")), p, width = 8, height = 5, dpi = 300)
  invisible(p)
}

## ============================================================================
## Helper: boxplot of a given cell type by subtype (Kruskal-Wallis + CS1 vs CS3)
## ============================================================================
plot_cell_box <- function(data_df, gene, title, cols, out_dir, labels,
                          group_levels = c("CS1", "CS2", "CS3", "CS4"),
                          subtitle = NULL) {
  fml  <- as.formula(paste0("`", gene, "` ~ group"))
  pval <- compare_means(fml, data = data_df, method = "kruskal.test")$p.format
  ymax <- max(data_df[[gene]], na.rm = TRUE)

  data_df$group <- factor(data_df$group, levels = group_levels)
  p <- ggplot(data_df, aes(x = group, y = .data[[gene]], fill = group, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5,
                 colour = "black", size = 0.4) +
    geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    theme_minimal(base_size = 15) +
    labs(x = "", y = "Cell Proportion", title = title, subtitle = subtitle) +
    theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.title.x    = element_text(size = 16, face = "bold"),
          axis.title.y    = element_text(size = 14),
          axis.text.x     = element_text(size = 14, face = "bold"),
          axis.text.y     = element_text(size = 14),
          legend.position = "none",
          plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle   = element_text(hjust = 0.5, size = 14, face = "bold")) +
    coord_cartesian(clip = "off", ylim = c(min(data_df[[gene]], na.rm = TRUE), ymax * 1.25)) +
    annotation_custom(
      grob = textGrob(label = paste0("p = ", pval), gp = gpar(fontsize = 15), hjust = 0.5),
      xmin = 1.1, xmax = 1.5, ymin = ymax * 0.85, ymax = ymax * 1.05) +
    stat_compare_means(aes(group = group),
                       method      = "wilcox.test",
                       label       = "p.format",
                       comparisons = list(c("CS1", "CS3")),
                       label.y     = ymax * 0.75)

  ggsave(file.path(out_dir, paste0(labels, gene, "_plot.png")), p, width = 4, height = 4, dpi = 300)
  ggsave(file.path(out_dir, paste0(labels, gene, "_plot.pdf")), p, width = 4, height = 4)
  invisible(p)
}

## ============================================================================
## Helper: per-cell-type boxplots loop (writes to <out_dir>/plot_test)
## ============================================================================
plot_all_celltypes <- function(ciber.res, group, out_dir, labels, cols, subtitle = NULL) {
  plot_dir <- file.path(out_dir, "plot_test")
  if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

  ann_df <- data.frame(cell = colnames(ciber.res),
                       labels = gsub("[^[:alnum:]]", "_", colnames(ciber.res)))
  rownames(ann_df) <- ann_df$labels

  data_df <- as.data.frame(ciber.res)
  data_df$group <- group$clust
  colnames(data_df) <- gsub("[^[:alnum:]]", "_", colnames(data_df))
  cell_list <- colnames(data_df)[1:(ncol(data_df) - 1)]

  for (gene in cell_list) {
    plot_cell_box(data_df, gene, title = ann_df[gene, "cell"], cols = cols,
                  out_dir = plot_dir, labels = labels, subtitle = subtitle)
  }
}

## ============================================================================
## Helper: per-subtype boxplots for ESTIMATE / MCPcounter scores (pairwise Wilcoxon)
## ============================================================================
plot_score_panel <- function(data_df, score_cols, cols, out_dir) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  group_levels <- unique(data_df$group)
  if (length(cols) < length(group_levels)) cols <- rep(cols, length.out = length(group_levels))
  names(cols) <- group_levels

  for (pathway in score_cols) {
    kw_p <- kruskal.test(reformulate("group", response = pathway), data = data_df)$p.value
    pairwise_res <- data_df |>
      pairwise_wilcox_test(formula = as.formula(paste(pathway, "~ group")),
                           p.adjust.method = "BH") |>
      filter(p.adj < 0.05) |>
      add_significance("p.adj") |>
      add_xy_position(x = "group", step.increase = 0.18)

    p <- ggplot(data_df, aes(x = group, y = .data[[pathway]], fill = group, colour = group)) +
      geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5,
                   colour = "black", linewidth = 0.4) +
      geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
      scale_fill_manual(values = cols) +
      scale_colour_manual(values = cols) +
      theme_minimal(base_size = 15) +
      labs(x = "", y = "", title = pathway,
           subtitle = paste0("Kruskal-Wallis p = ", signif(kw_p, 3))) +
      theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
            axis.title.x    = element_text(size = 16, face = "bold"),
            axis.title.y    = element_text(size = 14),
            axis.text.x     = element_text(size = 14, face = "bold"),
            axis.text.y     = element_text(size = 14),
            legend.position = "none",
            plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
            plot.subtitle   = element_text(hjust = 0, size = 12)) +
      coord_cartesian(clip = "off",
                      ylim = c(min(data_df[[pathway]], na.rm = TRUE),
                               max(data_df[[pathway]], na.rm = TRUE) * 1.6))

    if (nrow(pairwise_res) > 0) {
      p <- p + stat_pvalue_manual(pairwise_res, label = "p.adj.signif",
                                  tip.length = 0.01, hide.ns = TRUE, size = 5)
    }
    ggsave(file.path(out_dir, paste0(pathway, "_plot.png")), p, width = 4, height = 4, dpi = 300)
    ggsave(file.path(out_dir, paste0(pathway, "_plot.pdf")), p, width = 4, height = 4)
  }
}

## ############################################################################
## 1. Discovery cohort — CIBERSORT
## ############################################################################
load(file.path(ORIGINAL_DIR, "00_rawdata_2", "exp_merge_combat_nonneg.rData"))
exp_MASLD <- as.data.frame(res)

group <- read.csv(file.path(ORIGINAL_DIR, "01_MOVICS", "01_MOVICS_cluster_results.csv"), row.names = 1)
group$clust <- paste0("CS", group$clust)

ciber.res <- run_cibersort(exp_MASLD, output, tag = "11.")
ciber.res <- read.csv(file.path(output, paste0(labels, "11.CIBERSORT_ciber.res.csv")), row.names = 1)

plot_composition(ciber.res, group, output, labels, tag = "11.")

## ############################################################################
## 2. Discovery cohort — stacked boxplots (two rows, first/second 11 cell types)
## ############################################################################
ciber_res_long <- melt({ tmp <- as.data.frame(ciber.res); tmp$group <- group$clust; tmp },
                       id.vars = c("group"))
colnames(ciber_res_long) <- c("group", "cell type", "Relative percentage")
ciber_res_long$cell_type_num <- as.numeric(factor(ciber_res_long$`cell type`,
                                                  levels = unique(ciber_res_long$`cell type`)))
ciber_res_long$`cell type` <- str_replace_all(ciber_res_long$`cell type`, "\\.", " ")

make_half_panel <- function(df, cell_type_num_le, y_max, cols) {
  ggplot(df, aes(x = fct_reorder(`cell type`, `Relative percentage`),
                 y = `Relative percentage`, fill = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5,
                 colour = "black", size = 0.4) +
    geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = cols) +
    theme_minimal(base_size = 15) +
    labs(x = NULL, y = "Estimated Proportion") +
    theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.title.x    = element_text(size = 16, face = "bold"),
          axis.title.y    = element_text(size = 14),
          axis.text.x     = element_text(angle = 60, hjust = 1, size = 10, face = "bold"),
          axis.text.y     = element_text(size = 12),
          legend.position = "top",
          panel.grid.minor = element_blank()) +
    stat_compare_means(aes(group = group, label = ..p.signif..),
                       method = "kruskal.test", label.y = y_max * 0.9)
}

p1 <- make_half_panel(subset(ciber_res_long, cell_type_num <= 11),
                      y_max = max(subset(ciber_res_long, cell_type_num <= 11)$`Relative percentage`, na.rm = TRUE),
                      cols = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"))
p2 <- make_half_panel(subset(ciber_res_long, cell_type_num > 11),
                      y_max = max(subset(ciber_res_long, cell_type_num > 11)$`Relative percentage`, na.rm = TRUE),
                      cols = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"))
p3 <- ggarrange(p1, p2, nrow = 2, ncol = 1, common.legend = TRUE, legend = "top")
ggsave(file.path(output, paste0(labels, "result02.pdf")), p3, width = 8, height = 8)
ggsave(file.path(output, paste0(labels, "result02.png")), p3, width = 8, height = 8, dpi = 300)

## ############################################################################
## 3. Discovery cohort — per-cell-type boxplots
## ############################################################################
plot_all_celltypes(ciber.res, group, output, labels, cols = disc_cols)

## ############################################################################
## 4. External cohorts — CIBERSORT + composition + stacked boxes + per-cell plots
## ############################################################################
for (GSE in gse_list) {
  out_dir <- file.path(output, GSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  exp <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(GSE, "_exp.csv")), row.names = 1)
  grp <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(GSE, "_group.csv")), row.names = 1)
  cg  <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(labels, GSE, "cluster_res.csv")), row.names = 1)
  cg$clust <- paste0("CS", cg$clust)
  cg <- cg[rownames(grp), ]
  grp$clust <- cg$clust

  ciber.res <- run_cibersort(exp, out_dir, tag = paste0(GSE, ""))
  ciber.res <- read.csv(file.path(out_dir, paste0(labels, GSE, "CIBERSORT_ciber.res.csv")), row.names = 1)

  plot_composition(ciber.res, grp, out_dir, labels, tag = GSE, title = GSE)

  # Stacked boxplots (2 rows)
  ciber_res_long <- melt({ tmp <- as.data.frame(ciber.res); tmp$group <- grp$clust; tmp },
                         id.vars = c("group"))
  colnames(ciber_res_long) <- c("group", "cell type", "Relative percentage")
  ciber_res_long$cell_type_num <- as.numeric(factor(ciber_res_long$`cell type`,
                                                    levels = unique(ciber_res_long$`cell type`)))
  ciber_res_long$`cell type` <- str_replace_all(ciber_res_long$`cell type`, "\\.", " ")

  p1 <- make_half_panel(subset(ciber_res_long, cell_type_num <= 11),
                        y_max = max(subset(ciber_res_long, cell_type_num <= 11)$`Relative percentage`, na.rm = TRUE),
                        cols = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"))
  p2 <- make_half_panel(subset(ciber_res_long, cell_type_num > 11),
                        y_max = max(subset(ciber_res_long, cell_type_num > 11)$`Relative percentage`, na.rm = TRUE),
                        cols = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"))
  legend <- get_legend(
    ggplot(ciber_res_long, aes(x = `cell type`, y = `Relative percentage`, fill = group)) +
      geom_boxplot() +
      scale_fill_manual(values = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"), name = "Group") +
      theme(legend.position = "top",
            legend.title    = element_text(size = 14, face = "bold"),
            legend.text     = element_text(size = 12)))
  p3 <- ggarrange(p1, p2, nrow = 2, ncol = 1, common.legend = TRUE,
                  legend = "top", legend.grob = legend)
  p3 <- annotate_figure(p3, top = text_grob("Cell Type Distribution by Group",
                                            face = "bold", size = 18))
  ggsave(file.path(out_dir, paste0(labels, GSE, "result02.pdf")), p3, width = 8, height = 8)
  ggsave(file.path(out_dir, paste0(labels, GSE, "result02.png")), p3, width = 8, height = 8, dpi = 300)

  plot_all_celltypes(ciber.res, grp, out_dir, labels,
                     cols = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488"),
                     subtitle = GSE)
}

## ############################################################################
## 5. ESTIMATE scores (discovery)
## ############################################################################
exp_MASLD <- as.matrix(exp_MASLD)
storage.mode(exp_MASLD) <- "numeric"

# Drop all-NA / all-zero genes, then log2-transform
exp_MASLD <- exp_MASLD[rowSums(is.na(exp_MASLD)) < ncol(exp_MASLD), , drop = FALSE]
exp_MASLD <- exp_MASLD[rowSums(exp_MASLD, na.rm = TRUE) > 0, , drop = FALSE]
exp_MASLD <- log2(exp_MASLD + 1)

# Collapse duplicated gene symbols by mean
if (any(duplicated(rownames(exp_MASLD)))) {
  exp_MASLD <- as.data.frame(exp_MASLD) |>
    rownames_to_column("gene") |>
    group_by(gene) |>
    summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") |>
    column_to_rownames("gene") |>
    as.matrix()
}
rownames(exp_MASLD) <- trimws(rownames(exp_MASLD))

estimate_res <- deconvo_estimate(eset = exp_MASLD, project = "MASLD", platform = "affymetrix")
rownames(estimate_res) <- estimate_res$ID
estimate_res$group    <- group[rownames(estimate_res), "clust"]

estimate_dir <- file.path(output, "ESTIMATE")
if (!dir.exists(estimate_dir)) dir.create(estimate_dir, recursive = TRUE)
write.csv(estimate_res, file.path(estimate_dir, "MASLD_ESTIMATE_scores.csv"), row.names = FALSE)

plot_score_panel(estimate_res, score_cols = colnames(estimate_res)[3:6],
                 cols = disc_cols, out_dir = estimate_dir)

## ############################################################################
## 6. MCPcounter scores (discovery)
## ############################################################################
mcp_res <- deconvo_mcpcounter(eset = exp_MASLD, project = "MASLD")
rownames(mcp_res) <- estimate_res$ID
mcp_res$group <- group[rownames(mcp_res), "clust"]

mcp_dir <- file.path(output, "MCPcounter")
if (!dir.exists(mcp_dir)) dir.create(mcp_dir, recursive = TRUE)
write.csv(mcp_res, file.path(mcp_dir, "MCPcounter_res.csv"))

plot_score_panel(mcp_res, score_cols = "Fibroblasts_MCPcounter",
                 cols = disc_cols, out_dir = mcp_dir)
