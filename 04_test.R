rm(list = ls()); gc()

## -------- Paths & directories --------
ORIGINAL_DIR <- ""                        # working directory (set as needed)
output       <- file.path(ORIGINAL_DIR, "04_test")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
setwd(ORIGINAL_DIR)

# Purpose:
#   (1) Apply the 40-marker classifier to external validation cohorts (NTP/PAM),
#       report agreement, and plot clinical features by predicted subtype.
#   (2) Run PCoA on each validation cohort to assess separation of the 4 subtypes.
# All distances are Euclidean; PERMANOVA uses 999 permutations.

## -------- Libraries --------
library(MOVICS)
library(ggplot2)
library(dplyr)
library(ggpubr)
library(tidyr)
library(tidyverse)
library(caret)
library(ggplotify)
library(grid)
library(gtable)
library(vegan)
library(ape)
library(pairwiseAdonis)

## -------- Settings --------
labels       <- "04_test_"
dist_method  <- "euclidean"
permutations <- 999

cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
          "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

## -------- Shared inputs --------
load(file.path(ORIGINAL_DIR, "00_rawdata_2", "exp_merge_combat_nonneg.rData"))
exp_MASLD <- as.data.frame(res)
load(file.path(ORIGINAL_DIR, "00_KEGG_get", "ssgsea_01.rData"))

enzyme_list <- read.csv(file.path(ORIGINAL_DIR, "00_rawdata",
                                  "metabolism_subtyping_genes", "gendis_101589_mmc2.csv"))
enzyme_gene <- unique(enzyme_list$Genes)
enzyme_gene <- intersect(enzyme_gene, rownames(exp_MASLD))
mat <- as.matrix(exp_MASLD[enzyme_gene, ])

# 40-marker template from the discovery run
marker.up <- readRDS(file.path(ORIGINAL_DIR, "01_MOVICS", "marker.up.rData"))

# Consensus clustering result used by PAM
load(file.path(ORIGINAL_DIR, "01_MOVICS", "moic.res.list_02.rda"))
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = dist_method,
                               linkage       = "average")
save(cmoic.brca, file = file.path(output, "cmoic.brca.Rdata"))

## ============================================================================
## Helper 1: NTP + PAM prediction for one cohort, save inputs & cluster labels
## ============================================================================
predict_cohort <- function(GSE, exp, pdata, group, outdir) {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

  # Persist the cleaned inputs for downstream steps
  write.csv(exp,   file.path(outdir, paste0(GSE, "_exp.csv")))
  write.csv(group, file.path(outdir, paste0(GSE, "_group.csv")))
  write.csv(pdata, file.path(outdir, paste0(GSE, "_pdata.csv")))

  set.seed(123)
  ntp.pred <- runNTP(expr       = exp,
                     templates  = marker.up$templates,
                     scaleFlag  = TRUE,
                     centerFlag = TRUE,
                     doPlot     = TRUE,
                     fig.name   = paste0("NTP HEATMAP FOR ", GSE),
                     width      = 8,
                     height     = 6)

  pam.pred <- runPAM(train.expr = as.data.frame(mat),
                     moic.res   = cmoic.brca,
                     test.expr  = as.data.frame(exp))

  runKappa(subt1     = ntp.pred$clust.res$clust,
           subt2     = pam.pred$clust.res$clust,
           subt1.lab = "NTP",
           subt2.lab = "PAM",
           fig.name  = paste0("CONSISTENCY HEATMAP FOR ", GSE),
           fig.path  = outdir)

  cluster_res <- ntp.pred$clust.res
  write.csv(cluster_res, file.path(outdir, paste0(labels, GSE, "cluster_res.csv")))
  cluster_res
}

## ============================================================================
## Helper 2: boxplot of a continuous clinical variable by subtype (Kruskal-Wallis)
## ============================================================================
plot_box_by_cluster <- function(df, var, cluster_col = "Cluster",
                                cols, outdir, prefix, ymax_frac = 1.05) {
  fml  <- as.formula(paste(var, "~", cluster_col))
  pval <- compare_means(fml, data = df, method = "kruskal.test")$p.format
  ymax <- max(df[[var]], na.rm = TRUE)

  p <- ggplot(df, aes(x = .data[[cluster_col]], y = .data[[var]],
                      fill = .data[[cluster_col]], colour = .data[[cluster_col]])) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5,
                 colour = "black", size = 0.4) +
    geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    theme_minimal(base_size = 15) +
    labs(x = "", y = var, title = "") +
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
      ymin = ymax * (ymax_frac - 0.03), ymax = ymax * ymax_frac)

  ggsave(file.path(outdir, paste0(prefix, ".png")), p, width = 4, height = 4, dpi = 300)
  ggsave(file.path(outdir, paste0(prefix, ".pdf")), p, width = 4, height = 4)
  p
}

## ============================================================================
## Helper 3: stacked bar of a categorical variable by subtype
## ============================================================================
plot_bar_by_cluster <- function(df, var, cluster_col = "Cluster",
                                cols, outdir, prefix,
                                legend_title = var, legend_pos = "top",
                                legend_title_size = 10, legend_text_size = 10,
                                width = 4, height = 4) {
  df$tmp_var <- factor(df[[var]])
  prop <- df %>%
    group_by(.data[[cluster_col]], tmp_var) %>%
    summarise(Count = n(), .groups = "drop") %>%
    group_by(.data[[cluster_col]]) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()

  p <- ggplot(prop, aes(x = .data[[cluster_col]], y = Proportion, fill = tmp_var)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6,
             alpha = 0.8, colour = "black", linewidth = 0.4) +
    scale_fill_manual(values = cols, name = legend_title) +
    theme_minimal(base_size = 15) +
    labs(x = "", y = "Proportion") +
    theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.title.x    = element_text(size = 16, face = "bold"),
          axis.title.y    = element_text(size = 14),
          axis.text.x     = element_text(size = 14, face = "bold"),
          axis.text.y     = element_text(size = 14),
          legend.position = legend_pos,
          legend.title    = element_text(size = legend_title_size, face = "bold"),
          legend.text     = element_text(size = legend_text_size),
          plot.title      = element_text(hjust = 0.5, size = 16, face = "bold")) +
    coord_cartesian(clip = "off")

  ggsave(file.path(outdir, paste0(prefix, ".png")), p, width = width, height = height, dpi = 300)
  ggsave(file.path(outdir, paste0(prefix, ".pdf")), p, width = width, height = height)
  p
}

## ============================================================================
## Helper 4: PCoA + PERMANOVA + pairwiseAdonis + 2D plot + pairwise table
## ============================================================================
run_pcoa <- function(dist_mat, meta, out_dir, labels, table_h = 3) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  adon_res     <- adonis2(dist_mat ~ clust, data = meta, permutations = permutations)
  p_adonis     <- adon_res$`Pr(>F)`[1]

  bd           <- betadisper(dist_mat, meta$clust)
  bd_test      <- permutest(bd, permutations = permutations)
  p_betadisper <- bd_test$tab[1, "Pr(>F)"]

  pw <- pairwiseAdonis::pairwise.adonis(x = as.matrix(dist_mat), factors = meta$clust,
                                        sim.method = dist_method, p.adjust.m = "BH",
                                        perm = permutations)
  pw_out <- if (is.data.frame(pw)) pw else
    tryCatch(as.data.frame(pw), error = function(e) data.frame(result = I(list(pw))))
  write.csv(pw_out, file.path(out_dir, "pairwiseAdonis_results.csv"),
            row.names = FALSE, quote = FALSE)

  pcoa_res <- pcoa(dist_mat, correction = "none")
  eig_perc <- round(100 * pcoa_res$values$Relative_eig, 2)

  scores <- as.data.frame(pcoa_res$vectors[, 1:2])
  colnames(scores) <- c("PCoA1", "PCoA2")
  scores$Sample <- rownames(scores)
  scores <- merge(scores, meta[, c("Sample", "clust")], by = "Sample", sort = FALSE)

  p <- ggplot(scores, aes(x = PCoA1, y = PCoA2, color = clust, fill = clust)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(aes(group = clust), geom = "polygon", alpha = 0.2,
                 color = "black", level = 0.95) +
    scale_color_manual(values = cols, name = "Group") +
    scale_fill_manual(values = cols, name = "Group") +
    labs(title    = paste0("PCoA (", dist_method, " distance)"),
         subtitle = paste0("PERMANOVA p=", signif(p_adonis, 3),
                           "; betadisper p=", signif(p_betadisper, 3)),
         x = paste0("PCoA1 (", eig_perc[1], "%)"),
         y = paste0("PCoA2 (", eig_perc[2], "%)")) +
    theme_minimal(base_size = 15) +
    theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.title.x    = element_text(size = 16, face = "bold"),
          axis.title.y    = element_text(size = 16, face = "bold"),
          axis.text.x     = element_text(size = 14, face = "bold"),
          axis.text.y     = element_text(size = 14),
          legend.position = "right",
          legend.title    = element_text(size = 14, face = "bold"),
          legend.text     = element_text(size = 12),
          plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle   = element_text(hjust = 0.5, size = 12)) +
    coord_cartesian(clip = "off")
  ggsave(file.path(out_dir, paste0(labels, "PCoA_pairwiseAdonis.png")), p, width = 6, height = 5, dpi = 300)
  ggsave(file.path(out_dir, paste0(labels, "PCoA_pairwiseAdonis.pdf")), p, width = 6, height = 5)

  # Pairwise comparison table with three-line styling
  pw_out <- pw_out[, c(1, 5, 6, 7, 8)]
  tb <- tableGrob(
    pw_out, rows = NULL,
    theme = ttheme_default(
      core    = list(fg_params = list(cex = 0.8),
                     bg_params = list(fill = c("white", "white"), col = NA)),
      colhead = list(fg_params = list(cex = 0.8, fontface = "bold"),
                     bg_params = list(fill = "white", col = NA))))
  tb <- gtable_add_grob(
    tb,
    grobs = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                         y0 = unit(c(1, -1, 0), "lines"), y1 = unit(c(1, -1, 0), "lines"),
                         gp = gpar(lwd = 1.5)),
    t = 1, b = 1, l = 1, r = ncol(pw_out))
  tb <- gtable_add_grob(
    tb,
    grobs = segmentsGrob(x0 = unit(0, "npc"), x1 = unit(1, "npc"),
                         y0 = unit(0, "lines"), y1 = unit(0, "lines"),
                         gp = gpar(lwd = 2.5)),
    t = nrow(pw_out) + 1, b = nrow(pw_out) + 1, l = 1, r = ncol(pw_out))

  grid.newpage()
  png(file.path(out_dir, paste0(labels, "PCoA_pairwiseAdonis_table.png")),
      width = 6, height = table_h, res = 300, units = "in")
  grid.draw(tb); dev.off()
  grid.newpage()
  pdf(file.path(out_dir, paste0(labels, "PCoA_pairwiseAdonis_table.pdf")),
      width = 6, height = table_h)
  grid.draw(tb); dev.off()
}

## ============================================================================
## Helper 5: load a cohort's predicted labels and run the PCoA panel
## ============================================================================
pcoa_for_cohort <- function(GSE) {
  base_dir <- file.path(output, GSE)
  cluster_group <- read.csv(file.path(base_dir, paste0(labels, GSE, "cluster_res.csv")), row.names = 1)
  exp          <- read.csv(file.path(base_dir, paste0(GSE, "_exp.csv")), row.names = 1)
  cluster_group <- cluster_group[colnames(exp), ]
  cluster_group$clust <- paste0("CS", cluster_group$clust)

  expr_mat <- as.matrix(exp)
  mode(expr_mat) <- "numeric"

  meta <- cluster_group
  meta$Sample <- meta$samID
  if (!all(colnames(expr_mat) %in% meta$Sample)) stop("expr_mat columns must be contained in meta$Sample")
  meta <- meta[match(colnames(expr_mat), meta$Sample), , drop = FALSE]

  dist_mat <- dist(t(expr_mat), method = dist_method)
  run_pcoa(dist_mat, meta, base_dir, labels, table_h = 3)
}

## ############################################################################
## PART 1: External validation — NTP / PAM + clinical features
## ############################################################################

## -------- GSE193066 --------
GSE    <- "GSE193066"
outdir <- file.path(output, GSE)
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp_1.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)

exp <- exp[, -1]
annotation_map <- setNames(pdata$geo_accession, pdata$title)
colnames(exp)  <- ifelse(colnames(exp) %in% names(annotation_map),
                         annotation_map[colnames(exp)], colnames(exp))
group$group <- "NAFLD"

cluster_res <- predict_cohort(GSE, exp, pdata, group, outdir)
table(cluster_res$clust)

p_group <- pdata[, 42:47]
cluster_res <- cluster_res[rownames(p_group), ]
p_group$cluster <- paste0("CS", cluster_res$clust)
p_group <- p_group[, c("age.ch1", "biopsy.ch1", "fibrosis.stage.ch1",
                       "nas.ch1", "nafld.based.risk.prediction.ch1", "sex.ch1", "cluster")]
colnames(p_group) <- c("Age", "Biopsy", "Fibrosis_stage", "NAS",
                       "NAFLD_based_risk_prediction", "Sex", "Cluster")

plot_box_by_cluster(p_group, "Age", cols = cols, outdir = outdir,
                    prefix = paste0(labels, GSE, "_Age_cluster_"), ymax_frac = 1.05)
plot_bar_by_cluster(p_group, "Sex", cols = cols, outdir = outdir,
                    prefix = paste0(labels, GSE, "_Sex_cluster_"),
                    legend_title = "Gender", legend_pos = "top")
plot_box_by_cluster(p_group, "Fibrosis_stage", cols = cols, outdir = outdir,
                    prefix = paste0(labels, GSE, "_Fibrosis_stage_cluster_box_"))
plot_bar_by_cluster(p_group, "Fibrosis_stage", cols = cols, outdir = outdir,
                    prefix = paste0(labels, GSE, "_Fibrosis_stage_cluster_bar_"),
                    legend_title = "Fibrosis Stage", legend_pos = "top")
plot_box_by_cluster(p_group, "NAS", cols = cols, outdir = outdir,
                    prefix = paste0(labels, GSE, "_NAS_cluster_box_"), ymax_frac = 1.05)
plot_bar_by_cluster(p_group, "NAS", cols = cols, outdir = outdir,
                    prefix = paste0(labels, GSE, "_NAS_cluster_bar_"),
                    legend_title = "NAS", legend_pos = "top")

p_group2 <- subset(p_group, NAFLD_based_risk_prediction != "")
plot_bar_by_cluster(p_group2, "NAFLD_based_risk_prediction", cols = cols,
                    outdir = outdir,
                    prefix = paste0(labels, GSE, "_NAFLD_based_risk_prediction_cluster_bar_"),
                    legend_title = "Risk Prediction", legend_pos = "top")

## -------- GSE213621 --------
GSE    <- "GSE213621"
outdir <- file.path(output, GSE)
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
group  <- subset(group, group == "NAFLD")
exp    <- exp[, rownames(group)]
pdata  <- pdata[rownames(group), ]

cluster_res <- predict_cohort(GSE, exp, pdata, group, outdir)
p_group <- pdata[, 37, drop = FALSE]
cluster_res <- cluster_res[rownames(p_group), ]
p_group$cluster <- paste0("CS", cluster_res$clust)
colnames(p_group) <- c("Fibrosis_stage", "Cluster")

plot_bar_by_cluster(p_group, "Fibrosis_stage", cols = cols, outdir = outdir,
                    prefix = paste0(labels, GSE, "_Fibrosis_stage_cluster_bar_"),
                    legend_title = "Fibrosis Stage", legend_pos = "top")

## -------- GSE83452 --------
GSE    <- "GSE83452"
outdir <- file.path(output, GSE)
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
group  <- subset(group, group != "Obese")
exp    <- exp[, rownames(group)]
pdata  <- pdata[rownames(group), ]

cluster_res <- predict_cohort(GSE, exp, pdata, group, outdir)
p_group <- pdata[, 36:38, drop = FALSE]
cluster_res <- cluster_res[rownames(p_group), ]
p_group$cluster <- paste0("CS", cluster_res$clust)
colnames(p_group) <- c("Age", "Gender", "liver_status", "cluster")

p_group$liver_status <- ifelse(p_group$liver_status == "NASH", "MASH", "Non-MASH")

plot_box_by_cluster(p_group, "Age", cluster_col = "cluster", cols = cols,
                    outdir = outdir, prefix = paste0(labels, GSE, "_Age_cluster_"),
                    ymax_frac = 1.05)
plot_bar_by_cluster(p_group, "Gender", cluster_col = "cluster", cols = cols,
                    outdir = outdir, prefix = paste0(labels, GSE, "_Sex_cluster_"),
                    legend_title = "Gender", legend_pos = "right",
                    width = 5, height = 4)
plot_bar_by_cluster(p_group, "liver_status", cluster_col = "cluster", cols = cols,
                    outdir = outdir, prefix = paste0(labels, GSE, "_liver_status_cluster_"),
                    legend_title = "Liver Status", legend_pos = "right",
                    width = 5, height = 4)

## -------- GSE164760 --------
GSE    <- "GSE164760"
outdir <- file.path(output, GSE)
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
group  <- group %>%
  mutate(group = case_when(
    group %in% c("Healthy liver") ~ "Healthy",
    group %in% c("NASH liver") ~ "MASH",
    group %in% c("Non-tumoral NASH liver adjacent to HCC") ~ "AD-MASH",
    group %in% c("NASH-HCC tumor") ~ "NASH-HCC",
    group %in% c("Cirrhotic liver") ~ "Ci-liver",
    TRUE ~ NA_character_))
group <- subset(group, group != "Healthy")
exp   <- exp[, rownames(group)]
pdata <- pdata[rownames(group), ]

cluster_res <- predict_cohort(GSE, exp, pdata, group, outdir)
p_group <- pdata[, 29, drop = FALSE]
cluster_res <- cluster_res[rownames(p_group), ]
p_group$cluster <- paste0("CS", cluster_res$clust)
colnames(p_group) <- c("liver_tissue", "cluster")
p_group$liver_tissue <- group$group
p_group <- subset(p_group, liver_tissue != "NASH-HCC")

plot_bar_by_cluster(p_group, "liver_tissue", cluster_col = "cluster", cols = cols,
                    outdir = outdir, prefix = paste0(labels, GSE, "_liver_tissue_cluster_"),
                    legend_title = "Liver Tissue", legend_pos = "right",
                    width = 5, height = 4)
plot_bar_by_cluster(p_group, "cluster", cluster_col = "liver_tissue", cols = cols,
                    outdir = outdir, prefix = paste0(labels, GSE, "_cluster_liver_tissue_"),
                    legend_title = "Cluster", legend_pos = "right",
                    width = 5, height = 4)

## ############################################################################
## PART 2: PCoA on each validation cohort
## ############################################################################
for (GSE in c("GSE193066", "GSE213621", "GSE83452")) {
  pcoa_for_cohort(GSE)
}
