rm(list = ls()); gc()

## -------- Paths & directories --------
ORIGINAL_DIR <- ""                        # working directory (set as needed)
output       <- file.path(ORIGINAL_DIR, "03_PCoA")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
setwd(ORIGINAL_DIR)

# Purpose: PCoA ordination to visually assess separation of MASLD subtypes (k = 4)
# and to compare separation across data types (expression, metabolic, Fenton,
# stress) and clinical groups.

library(vegan)
library(ape)
library(ggplot2)
library(ggsci)
library(ggrepel)
library(pairwiseAdonis)
library(gridExtra)
library(grid)
library(gtable)
library(tidyverse)

dist_method  <- "euclidean"
permutations <- 999

## -------- Shared palettes --------
cols_clust <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
                "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")
cols_group <- c("#F39B7F", "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

## -------- Load expression matrix & cluster assignments --------
load(file.path(ORIGINAL_DIR, "00_rawdata_2", "exp_merge_combat_nonneg.rData"))
exp_MASLD <- as.data.frame(res)

group  <- read.csv(file.path(ORIGINAL_DIR, "01_MOVICS", "01_MOVICS_cluster_results.csv"), row.names = 1)
expr_mat <- as.matrix(exp_MASLD)
mode(expr_mat) <- "numeric"

meta <- group
meta$Sample <- meta$samID
meta$clust  <- as.character(meta$clust)
meta$clust  <- ifelse(meta$clust %in% c("1", "2", "3", "4"), paste0("CS", meta$clust), meta$clust)
table(meta$clust)
if (!all(colnames(expr_mat) %in% meta$Sample)) stop("expr_mat columns must be contained in meta$Sample")
meta <- meta[match(colnames(expr_mat), meta$Sample), , drop = FALSE]

dist_mat <- dist(t(expr_mat), method = dist_method)

## -------- Helper: PCoA + PERMANOVA + pairwiseAdonis + 2D/3D plots + table --------
run_pcoa <- function(dist_mat, meta, group_col, cols, out_dir, labels,
                     perm = permutations, dist_method = "euclidean",
                     table_h = 2.5) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  ## PERMANOVA + dispersion homogeneity
  adon_res     <- adonis2(dist_mat ~ meta[[group_col]], data = meta, permutations = perm)
  p_adonis     <- adon_res$`Pr(>F)`[1]

  bd           <- betadisper(dist_mat, meta[[group_col]])
  bd_test      <- permutest(bd, permutations = perm)
  p_betadisper <- bd_test$tab[1, "Pr(>F)"]

  ## Pairwise PERMANOVA with BH correction
  pw <- pairwiseAdonis::pairwise.adonis(x = as.matrix(dist_mat), factors = meta[[group_col]],
                                        sim.method = dist_method, p.adjust.m = "BH", perm = perm)
  pw_out <- if (is.data.frame(pw)) pw else
    tryCatch(as.data.frame(pw), error = function(e) data.frame(result = I(list(pw))))
  write.csv(pw_out, file = file.path(out_dir, "pairwiseAdonis_results.csv"),
            row.names = FALSE, quote = FALSE)

  ## PCoA
  pcoa_res <- pcoa(dist_mat, correction = "none")
  eig_perc <- round(100 * pcoa_res$values$Relative_eig, 2)

  scores <- as.data.frame(pcoa_res$vectors[, 1:2])
  colnames(scores) <- c("PCoA1", "PCoA2")
  scores$Sample <- rownames(scores)
  scores <- merge(scores, meta[, c("Sample", group_col)], by = "Sample", sort = FALSE)
  colnames(scores)[ncol(scores)] <- "grp"

  ## 2D plot
  p <- ggplot(scores, aes(x = PCoA1, y = PCoA2, color = grp, fill = grp)) +
    geom_point(size = 2.5, alpha = 0.8) +
    stat_ellipse(aes(group = grp), geom = "polygon", alpha = 0.2, color = "black", level = 0.95) +
    scale_color_manual(values = cols, name = "Group") +
    scale_fill_manual(values = cols, name = "Group") +
    labs(
      title    = paste0("PCoA (", dist_method, " distance)"),
      subtitle = paste0("PERMANOVA p=", signif(p_adonis, 3),
                        "; betadisper p=", signif(p_betadisper, 3)),
      x = paste0("PCoA1 (", eig_perc[1], "%)"),
      y = paste0("PCoA2 (", eig_perc[2], "%)")
    ) +
    theme_minimal(base_size = 15) +
    theme(
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.title.x    = element_text(size = 16, face = "bold"),
      axis.title.y    = element_text(size = 16, face = "bold"),
      axis.text.x     = element_text(size = 14, face = "bold"),
      axis.text.y     = element_text(size = 14),
      legend.position = "right",
      legend.title    = element_text(size = 14, face = "bold"),
      legend.text     = element_text(size = 12),
      plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
      plot.subtitle   = element_text(hjust = 0.5, size = 12)
    ) +
    coord_cartesian(clip = "off")
  ggsave(file.path(out_dir, paste0(labels, "PCoA_cluster.png")), p, width = 6, height = 5, dpi = 300)
  ggsave(file.path(out_dir, paste0(labels, "PCoA_cluster.pdf")), p, width = 6, height = 5)

  ## 3D interactive plot
  scores3d <- as.data.frame(pcoa_res$vectors[, 1:3])
  colnames(scores3d) <- c("PCoA1", "PCoA2", "PCoA3")
  scores3d$Sample <- rownames(scores3d)
  scores3d <- merge(scores3d, meta[, c("Sample", group_col)], by = "Sample", sort = FALSE)
  colnames(scores3d)[ncol(scores3d)] <- "grp"
  eig_perc3 <- round(100 * pcoa_res$values$Relative_eig[1:3], 2)
  cols_use  <- setNames(cols[seq_along(unique(scores3d$grp))], unique(scores3d$grp))

  p3d <- plot_ly(
    data = scores3d, x = ~PCoA1, y = ~PCoA2, z = ~PCoA3,
    color = ~grp, colors = cols_use,
    type = "scatter3d", mode = "markers",
    marker = list(size = 5, opacity = 0.85),
    text = ~paste0("Sample: ", Sample, "<br>Group: ", grp,
                   "<br>PCoA1: ", round(PCoA1, 3),
                   "<br>PCoA2: ", round(PCoA2, 3),
                   "<br>PCoA3: ", round(PCoA3, 3)),
    hoverinfo = "text"
  ) %>%
    layout(
      title  = list(text = paste0("3D PCoA (", dist_method, " distance)"), x = 0.5),
      scene  = list(
        xaxis = list(title = paste0("PCoA1 (", eig_perc3[1], "%)")),
        yaxis = list(title = paste0("PCoA2 (", eig_perc3[2], "%)")),
        zaxis = list(title = paste0("PCoA3 (", eig_perc3[3], "%)"))
      ),
      legend = list(title = list(text = "Group"))
    )
  htmlwidgets::saveWidget(p3d, file.path(out_dir, paste0(labels, "PCoA_cluster_3D.html")),
                          selfcontained = TRUE)

  ## Pairwise table
  pw_out <- pw_out[, c(1, 5, 6, 7, 8)]
  make_table <- function(pw_out, bottom_row) {
    tb <- tableGrob(
      pw_out, rows = NULL,
      theme = ttheme_default(
        core    = list(fg_params = list(cex = 0.8),
                       bg_params = list(fill = c("white", "white"), col = NA)),
        colhead = list(fg_params = list(cex = 0.8, fontface = "bold"),
                       bg_params = list(fill = "white", col = NA))
      ))
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
      t = bottom_row, b = bottom_row, l = 1, r = ncol(pw_out))
    tb
  }
  table <- make_table(pw_out, nrow(pw_out) + 1)

  grid.newpage()
  png(file.path(out_dir, paste0(labels, "PCoA_cluster_table.png")),
      width = 6, height = table_h, res = 300, units = "in")
  grid.draw(table); dev.off()
  grid.newpage()
  pdf(file.path(out_dir, paste0(labels, "PCoA_cluster_table.pdf")),
      width = 6, height = table_h)
  grid.draw(table); dev.off()
}

## ============================================================================
## 1. Subtype separation on full transcriptome
## ============================================================================
run_pcoa(dist_mat, meta, "clust", cols_clust,
         out_dir = file.path(output, "pcoa_pairwiseAdonis_cluster"),
         labels  = "07_PCoA_", table_h = 2.5)

## ============================================================================
## 2. Subtype separation vs. clinical group (SS / MASH)
## ============================================================================
group_raw <- read.csv(file.path(ORIGINAL_DIR, "00_rawdata_2", "00_rawdata_group_merged.csv"), row.names = 1)
group_MASLD <- group_raw[colnames(exp_MASLD), ]
table(group_MASLD$group)
group_MASLD <- group_MASLD %>%
  mutate(group = case_when(
    group %in% c("MASLD", "NAFLD") ~ "MASLD",
    group %in% c("NAFL") ~ "SS",
    group %in% c("NASH") ~ "MASH",
    TRUE ~ NA_character_))
table(group_MASLD$group)
group_MASLD <- group_MASLD[rownames(meta), ]
group_MASLD$cluster <- as.character(meta$clust)
meta$group <- group_MASLD$group

run_pcoa(dist_mat, meta, "group", cols_group,
         out_dir = file.path(output, "pcoa_pairwiseAdonis_group"),
         labels  = "07_PCoA_", table_h = 1.5)

## ============================================================================
## 3. Separation based on metabolic pathway scores
## ============================================================================
load(file.path(ORIGINAL_DIR, "11_GSVA_k4", "11_GSVA_k4_gsva_results.Rdata"))
me_data <- as.data.frame(t(gsva_results))

meta_met <- meta[rownames(me_data), ]
i_list <- c("GOBP_AMINO_ACID_CATABOLIC_PROCESS", "GOBP_FATTY_ACID_CATABOLIC_PROCESS",
            "GOBP_GOLGI_VESICLE_TRANSPORT", "GOBP_DNA_CATABOLIC_PROCESS",
            "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY",
            "GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION",
            "GOBP_CARBOHYDRATE_CATABOLIC_PROCESS")
dist_met <- dist(me_data[, i_list], method = dist_method)

run_pcoa(dist_met, meta_met, "clust", cols_clust,
         out_dir = file.path(output, "metabolic_test"),
         labels  = "07_PCoA_", table_h = 2.5)

## ============================================================================
## 4. Separation based on compartment-resolved Fenton reaction scores
## ============================================================================
load(file.path(ORIGINAL_DIR, "16_Fenton", "16_Fenton_fenton_gsva_results.Rdata"))
me_data <- as.data.frame(t(gsva_results))
meta_fen <- meta[rownames(me_data), ]
dist_fen <- dist(me_data, method = dist_method)

run_pcoa(dist_fen, meta_fen, "clust", cols_clust,
         out_dir = file.path(output, "Fenton_test"),
         labels  = "07_PCoA_", table_h = 2.5)

## ============================================================================
## 5. Separation based on cellular stress scores
## ============================================================================
load(file.path(ORIGINAL_DIR, "13_stress_GSVA", "13_stress_GSVA_stress_gsva_results.Rdata"))
me_data <- as.data.frame(t(gsva_results))
meta_str <- meta[rownames(me_data), ]
dist_str <- dist(me_data, method = dist_method)

run_pcoa(dist_str, meta_str, "clust", cols_clust,
         out_dir = file.path(output, "stress_test"),
         labels  = "07_PCoA_", table_h = 2.5)
