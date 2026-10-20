rm(list = ls()); gc()

## -------- Paths & directories --------
ORIGINAL_DIR <- ""                        # working directory (set as needed)
output       <- file.path(ORIGINAL_DIR, "02_function")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
setwd(ORIGINAL_DIR)

# Purpose: quantify predefined functional categories across subtypes — build
# category-to-pathway maps, aggregate GSVA scores per sample, test one-vs-rest
# differences, and visualize as radar plots + ComplexHeatmap. k = 4.

## -------- Libraries --------
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)

## -------- Settings --------
labels       <- "02_function_"
group_levels <- c("CS1", "CS2", "CS3", "CS4")
group_colors <- c(CS1 = "#E41A1C", CS2 = "#377EB8", CS3 = "#4DAF4A", CS4 = "#FF7F00")

## ============================================================================
## Helper: load subtype labels as a sample-level data.frame
## ============================================================================
load_groups <- function(path) {
  group <- read.csv(path, row.names = 1)
  group$clust <- paste0("CS", group$clust)
  data.frame(sample = rownames(group), group = group$clust, stringsAsFactors = FALSE)
}

## ============================================================================
## Helper: one-vs-rest Wilcoxon test for each function × subtype
## ============================================================================
one_vs_rest_test <- function(gsva_cat_sample, func_order, group_levels) {
  bind_rows(lapply(func_order, function(fn) {
    df_f <- gsva_cat_sample %>%
      filter(Function == fn) %>%
      mutate(group = factor(group, levels = group_levels))

    bind_rows(lapply(group_levels, function(g) {
      x <- df_f %>% filter(group == g) %>% pull(Function_GSVA)
      y <- df_f %>% filter(group != g) %>% pull(Function_GSVA)
      if (length(x) < 2 || length(y) < 2) return(NULL)

      p <- tryCatch(stats::wilcox.test(x, y, exact = FALSE)$p.value,
                    error = function(e) NA_real_)
      tibble(Function = fn, group = g, p_value = p,
             median_group = median(x, na.rm = TRUE),
             median_other = median(y, na.rm = TRUE),
             direction    = ifelse(median(x, na.rm = TRUE) > median(y, na.rm = TRUE),
                                   "Higher", "Lower"))
    }))
  })) %>%
    mutate(fdr  = p.adjust(p_value, method = "BH"),
           star = case_when(fdr < 0.001 ~ "***",
                            fdr < 0.01  ~ "**",
                            fdr < 0.05  ~ "*",
                            TRUE ~ ""))
}

## ============================================================================
## Helper: safe Kruskal-Wallis p-value
## ============================================================================
safe_kw_p <- function(x, g) {
  tryCatch({
    if (length(unique(g[!is.na(x)])) < 2) return(NA_real_)
    if (length(unique(x[!is.na(x)])) < 2) return(NA_real_)
    stats::kruskal.test(x ~ g)$p.value
  }, error = function(e) NA_real_)
}

## ============================================================================
## Helper: radar plot for one subtype
## ============================================================================
make_radar_plot <- function(group_name, plot_df, pairwise_res, baseline_df,
                            y_lim, output_dir, function_map, color_map, stress_order) {
  n_axis <- length(stress_order)

  df_g <- plot_df %>%
    filter(group == group_name) %>%
    left_join(function_map, by = "Function") %>%
    select(Function_ID = ID, mean_GSVA) %>%
    complete(Function_ID = 1:n_axis, fill = list(mean_GSVA = 0)) %>%
    mutate(x = Function_ID) %>%
    arrange(Function_ID)
  df_g$Color <- color_map$Color[match(df_g$Function_ID, color_map$ID)]

  df_closed <- bind_rows(df_g, mutate(df_g[1, ], x = nrow(df_g) + 1))

  baseline_plot_df <- baseline_df %>%
    left_join(function_map, by = "Function") %>%
    select(Function_ID = ID, baseline_GSVA) %>%
    complete(Function_ID = 1:n_axis, fill = list(baseline_GSVA = 0)) %>%
    mutate(x = Function_ID) %>%
    arrange(Function_ID)
  baseline_closed <- bind_rows(baseline_plot_df,
                               mutate(baseline_plot_df[1, ], x = nrow(baseline_plot_df) + 1))

  sig_df <- pairwise_res %>%
    filter(group == group_name, star != "") %>%
    left_join(function_map, by = "Function") %>%
    distinct(ID, .keep_all = TRUE) %>%
    mutate(x = as.numeric(ID), y = y_lim[2] + 0.06 * diff(y_lim)) %>%
    rename(Function_ID = ID)

  p <- ggplot() +
    geom_path(data = baseline_closed, aes(x = x, y = baseline_GSVA, group = 1),
              color = "grey50", linetype = "dashed", linewidth = 1.0) +
    geom_path(data = df_closed, aes(x = x, y = mean_GSVA, group = 1),
              color = "#1A237E", linewidth = 1.5) +
    geom_point(data = df_g, aes(x = x, y = mean_GSVA, color = I(Color)),
               size = 5, stroke = 1, shape = 19) +
    geom_text(data = sig_df, aes(x = x, y = y, label = star),
              inherit.aes = FALSE, color = "red", size = 7, fontface = "bold") +
    coord_polar(start = -pi / 2, clip = "off") +
    scale_x_continuous(breaks = 1:n_axis, labels = 1:n_axis, limits = c(1, n_axis + 1)) +
    scale_y_continuous(limits = c(y_lim[1], y_lim[2] + 0.12 * diff(y_lim)),
                       breaks = pretty(y_lim, n = 5)) +
    labs(title = paste0(group_name, " Cell Stress GSVA Radar Plot"),
         x = NULL, y = "Raw GSVA") +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_line(color = "grey85", size = 0.5),
          axis.text.x  = element_text(size = 12, color = "black", face = "bold"),
          axis.text.y  = element_text(size = 10, color = "grey30"),
          axis.title.y = element_text(size = 12, face = "bold"),
          plot.title   = element_text(hjust = 0.5, face = "bold", size = 15),
          plot.margin  = margin(15, 20, 15, 20))

  ggsave(file.path(output_dir, paste0(group_name, "_cell_stress_radar.png")),
         p, width = 6, height = 6, dpi = 300)
  ggsave(file.path(output_dir, paste0(group_name, "_cell_stress_radar.pdf")),
         p, width = 6, height = 6)
  list(plot = p, data = df_g)
}

## ============================================================================
## Helper: ComplexHeatmap of functional GSVA scores with KW stars
## ============================================================================
draw_function_heatmap <- function(gsva_cat_sample, out_dir, group_levels,
                                  group_colors, file_prefix,
                                  title = "Functional GSVA Heatmap by CS Subtype") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  ann_df <- gsva_cat_sample %>%
    select(Sample, group) %>%
    distinct() %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    arrange(group, Sample) %>%
    column_to_rownames("Sample")
  sample_order <- rownames(ann_df)

  heatmap_df <- gsva_cat_sample %>%
    select(Function, Sample, Function_GSVA) %>%
    group_by(Function, Sample) %>%
    summarise(Function_GSVA = mean(Function_GSVA, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Sample, values_from = Function_GSVA) %>%
    column_to_rownames("Function")
  heatmap_mat <- as.matrix(heatmap_df[, sample_order, drop = FALSE])

  kw_results <- gsva_cat_sample %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    group_by(Function) %>%
    summarise(kw_p = safe_kw_p(Function_GSVA, group), .groups = "drop") %>%
    mutate(kw_fdr = p.adjust(kw_p, method = "BH"))
  sig_fun  <- kw_results %>% filter(kw_fdr < 0.05) %>% pull(Function)
  row_labels <- rownames(heatmap_mat)
  row_labels <- ifelse(row_labels %in% sig_fun, paste0(row_labels, "*"), row_labels)

  mat_min <- min(heatmap_mat, na.rm = TRUE)
  mat_max <- max(heatmap_mat, na.rm = TRUE)
  col_fun <- circlize::colorRamp2(c(mat_min, 0, mat_max),
                                  c("#2166AC", "white", "#B2182B"))

  ha <- HeatmapAnnotation(group = ann_df$group,
                          col = list(group = group_colors),
                          annotation_name_side = "left")

  ht <- Heatmap(
    heatmap_mat,
    name = "GSVA",
    col = col_fun,
    top_annotation = ha,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    column_order = sample_order,
    column_split = ann_df$group,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_labels = row_labels,
    row_names_gp = gpar(fontsize = 8),
    row_title = NULL,
    column_title = title,
    heatmap_legend_param = list(title = "Mean GSVA"),
    border = FALSE,
    use_raster = TRUE)

  pdf(file.path(out_dir, paste0(file_prefix, ".pdf")), width = 14, height = 10)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  png(file.path(out_dir, paste0(file_prefix, ".png")), width = 4200, height = 3000, res = 300)
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()

  write.csv(kw_results, file.path(out_dir, paste0(file_prefix, "_KW_results.csv")),
            row.names = FALSE)
}

## ############################################################################
## 1. Load inputs
## ############################################################################
pathway_set <- read.csv(file.path(ORIGINAL_DIR, "02_function", "pathways.csv"),
                        stringsAsFactors = FALSE, check.names = FALSE)$trems
pathway_set <- unique(na.omit(pathway_set))

group_df <- load_groups(file.path(ORIGINAL_DIR, "01_MOVICS", "01_MOVICS_cluster_results.csv"))

load(file.path(ORIGINAL_DIR, "02_GSVA_k4", "02_GSVA_k4_gsva_results.Rdata"))

## ############################################################################
## 2. Functional dictionary (category -> pathway list)
## ############################################################################
class_dict <- list(
  "Lipid biosynthesis" = c(
    "GOBP_LIPID_BIOSYNTHETIC_PROCESS", "Fatty acid synthesis",
    "Phospholipid synthesis-PA", "Phospholipid synthesis-PC", "Phospholipid synthesis-PE",
    "Phospholipid synthesis-PI", "Phospholipid synthesis-PS",
    "Triglyceride synthesis", "Ceramide synthesis"),
  "Fatty acid oxidation" = c(
    "GOBP_FATTY_ACID_BETA_OXIDATION", "GOBP_FATTY_ACID_CATABOLIC_PROCESS",
    "GOBP_FATTY_ACID_METABOLIC_PROCESS", "Beta-oxidation", "Triglyceride degradation"),
  "Cholesterol metabolism" = c(
    "GOBP_CHOLESTEROL_CATABOLIC_PROCESS",
    "GOBP_REGULATION_OF_CHOLESTEROL_BIOSYNTHETIC_PROCESS",
    "GOBP_REGULATION_OF_CHOLESTEROL_METABOLIC_PROCESS", "Mevalonate pathway"),
  "Bile acid metabolism" = c(
    "GOBP_BILE_ACID_BIOSYNTHETIC_PROCESS", "GOBP_BILE_ACID_METABOLIC_PROCESS"),
  "Mitochondrial function" = c(
    "GOBP_MITOCHONDRION_ORGANIZATION", "GOBP_AUTOPHAGY_OF_MITOCHONDRION"),
  "Oxidative phosphorylation" = c(
    "GOBP_OXIDATIVE_PHOSPHORYLATION", "GOBP_ATP_SYNTHESIS_COUPLED_ELECTRON_TRANSPORT",
    "GOBP_ATP_BIOSYNTHETIC_PROCESS", "GOBP_ATP_TRANSPORT"),
  "ROS" = c(
    "GOBP_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS", "GOBP_RESPONSE_TO_OXIDATIVE_STRESS",
    "GOBP_REGULATION_OF_CELLULAR_RESPONSE_TO_OXIDATIVE_STRESS",
    "GOBP_REGULATION_OF_RESPONSE_TO_OXIDATIVE_STRESS"),
  "Metabolic reprogramming" = c(
    "Arginine transportation", "Choline production", "Glutaminolysis",
    "Gluconeogenesis-specific", "Lysine degradation", "Proline synthesis",
    "Serine synthesis", "Tryptophan degradation", "Retinol biosynthesis",
    "Retinol metabolism", "Circadian rhythm", "N-glycosylation complex synthesis phase",
    "N-glycosylation initial phase", "N-glycosylation processing phase", "O-glycosylation",
    "Phospholipid degradation", "Reprogrammed lipid metabolism",
    "Pyrimidine RN de novo synthesis", "Pyrimidine RN degradation",
    "Pyrimidine RN salvage synthesis", "Purine RN de novo synthesis",
    "Purine RN degradation", "Purine RN salvage synthesis",
    "Purine dRN de novo synthesis", "Purine dRN salvage synthesis",
    "Pyrimidine dRN de novo synthesis", "Pyrimidine dRN salvage synthesis",
    "Fatty acid transporter", "Sialic acid synthesis", "Heparan sulfate synthesis",
    "Chondroitin sulfate synthesis", "Hyaluronic acid synthesis", "Keratan sulfate synthesis"),
  "Inflammatory response" = c("GOBP_INFLAMMATORY_RESPONSE"),
  "TNF signaling" = c(
    "GOBP_INFLAMMATORY_RESPONSE", "GOBP_CANONICAL_NF_KAPPAB_SIGNAL_TRANSDUCTION",
    "GOBP_REGULATION_OF_CHEMOKINE_MEDIATED_SIGNALING_PATHWAY"),
  "NF-kappaB signaling" = c("GOBP_CANONICAL_NF_KAPPAB_SIGNAL_TRANSDUCTION"),
  "Interferon signaling" = c(
    "GOBP_INTERFERON_ALPHA_PRODUCTION", "GOBP_INTERFERON_BETA_PRODUCTION",
    "GOBP_INTERFERON_MEDIATED_SIGNALING_PATHWAY",
    "GOBP_CELLULAR_RESPONSE_TO_INTERFERON_ALPHA",
    "GOBP_CELLULAR_RESPONSE_TO_INTERFERON_BETA",
    "GOBP_CELLULAR_RESPONSE_TO_TYPE_II_INTERFERON",
    "GOBP_CELLULAR_RESPONSE_TO_TYPE_III_INTERFERON",
    "GOBP_RESPONSE_TO_INTERFERON_ALPHA", "GOBP_RESPONSE_TO_INTERFERON_BETA",
    "GOBP_RESPONSE_TO_TYPE_I_INTERFERON", "GOBP_RESPONSE_TO_TYPE_II_INTERFERON",
    "GOBP_RESPONSE_TO_TYPE_III_INTERFERON",
    "GOBP_REGULATION_OF_TYPE_I_INTERFERON_MEDIATED_SIGNALING_PATHWAY",
    "GOBP_TYPE_I_INTERFERON_PRODUCTION", "GOBP_TYPE_II_INTERFERON_MEDIATED_SIGNALING_PATHWAY",
    "GOBP_TYPE_II_INTERFERON_PRODUCTION", "GOBP_TYPE_III_INTERFERON_PRODUCTION"),
  "ER stress" = c("GOBP_CELLULAR_RESPONSE_TO_UNFOLDED_PROTEIN"),
  "UPR"       = c("GOBP_CELLULAR_RESPONSE_TO_UNFOLDED_PROTEIN"),
  "Autophagy" = c(
    "GOBP_MACROAUTOPHAGY", "GOBP_AUTOPHAGY_OF_MITOCHONDRION", "GOBP_REGULATION_OF_AUTOPHAGY"),
  "Apoptosis" = c(
    "GOBP_EXECUTION_PHASE_OF_APOPTOSIS", "GOBP_REGULATION_OF_EXECUTION_PHASE_OF_APOPTOSIS"),
  "Extracellular matrix remodeling" = c(
    "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY", "GOBP_EXTRACELLULAR_MATRIX_CELL_SIGNALING",
    "GOBP_EXTRACELLULAR_TRANSPORT", "GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ASSEMBLY",
    "GOBP_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION",
    "GOBP_POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION"),
  "Fibrosis" = c(
    "GOBP_COLLAGEN_BIOSYNTHETIC_PROCESS", "GOBP_COLLAGEN_CATABOLIC_PROCESS",
    "GOBP_COLLAGEN_FIBRIL_ORGANIZATION", "GOBP_COLLAGEN_METABOLIC_PROCESS",
    "GOBP_REGULATION_OF_COLLAGEN_FIBRIL_ORGANIZATION",
    "GOBP_REGULATION_OF_COLLAGEN_METABOLIC_PROCESS"),
  "TGF-beta signaling" = c(
    "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA",
    "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION",
    "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA1_PRODUCTION",
    "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA2_PRODUCTION"),
  "EMT" = c(
    "GOBP_TRANSFORMING_GROWTH_FACTOR_BETA_RECEPTOR_SIGNALING_PATHWAY",
    "GOBP_RESPONSE_TO_TRANSFORMING_GROWTH_FACTOR_BETA",
    "GOBP_EXTRACELLULAR_MATRIX_ASSEMBLY", "GOBP_COLLAGEN_FIBRIL_ORGANIZATION",
    "GOBP_POSITIVE_REGULATION_OF_EXTRACELLULAR_MATRIX_ORGANIZATION"),
  "Immune cell infiltration" = c(
    "GOBP_LEUKOCYTE_AGGREGATION", "GOBP_LEUKOCYTE_CHEMOTAXIS", "GOBP_LEUKOCYTE_MIGRATION",
    "GOBP_REGULATION_OF_LEUKOCYTE_CHEMOTAXIS", "GOBP_REGULATION_OF_LEUKOCYTE_MIGRATION"),
  "Chemokines / Chemotaxis" = c(
    "GOBP_RESPONSE_TO_CHEMOKINE", "GOBP_REGULATION_OF_CHEMOKINE_MEDIATED_SIGNALING_PATHWAY"),
  "Antigen presentation" = c(
    "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
    "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN",
    "GOBP_REGULATION_OF_ANTIGEN_PROCESSING_AND_PRESENTATION")
)
func_order <- names(class_dict)

# Flat table: function -> pathway
df <- bind_rows(lapply(seq_along(class_dict), function(i) {
  data.frame(functions = names(class_dict)[i], pathway = class_dict[[i]])
}))
write.csv(df, file.path(output, "df_result_pathway.csv"))

# Annotate pathways with their gene members from the local c5 collection
load(file.path(output, "c5_geneset.Rdata"))
df2 <- data.frame(pathway = names(c5_sets),
                  genes   = sapply(c5_sets, function(x) paste(x, collapse = ", ")))
pathways_keep <- intersect(rownames(df2), df$pathway)
df3 <- df2[pathways_keep, ]
rownames(df) <- df$pathway
df_filtered <- df[df[, 2] %in% rownames(df2), ]
df_result <- merge(df_filtered, df3[, c("pathway", "genes")], by = "pathway", all.x = TRUE)
write.csv(df_result, file.path(output, "df_result_pathway_genes.csv"))

## ############################################################################
## 3. Aggregate GSVA to functional categories
## ############################################################################
gsva_df   <- as.data.frame(gsva_results) %>% tibble::rownames_to_column("Pathway")
gsva_long <- gsva_df %>%
  pivot_longer(cols = -Pathway, names_to = "Sample", values_to = "GSVA") %>%
  left_join(group_df, by = c("Sample" = "sample"))
if (any(is.na(gsva_long$group))) stop("Some samples do not match the group file.")

category_map <- bind_rows(lapply(func_order, function(fn) {
  pw <- intersect(class_dict[[fn]], gsva_df$Pathway)
  if (!length(pw)) return(NULL)
  data.frame(Pathway = pw, Function = fn, stringsAsFactors = FALSE)
})) %>% distinct(Pathway, Function)

gsva_cat_sample <- gsva_long %>%
  inner_join(category_map, by = "Pathway") %>%
  group_by(Sample, group, Function) %>%
  summarise(Function_GSVA = mean(GSVA, na.rm = TRUE), .groups = "drop")

gsva_cat_group <- gsva_cat_sample %>%
  group_by(group, Function) %>%
  summarise(mean_GSVA = mean(Function_GSVA, na.rm = TRUE), .groups = "drop")

baseline_df <- gsva_cat_sample %>%
  group_by(Function) %>%
  summarise(baseline_GSVA = mean(Function_GSVA, na.rm = TRUE), .groups = "drop")

## ############################################################################
## 4. One-vs-rest significance test
## ############################################################################
pairwise_res <- one_vs_rest_test(gsva_cat_sample, func_order, group_levels)

output_radar <- file.path(output, "Radar")
if (!dir.exists(output_radar)) dir.create(output_radar, recursive = TRUE)
write.csv(pairwise_res,     file.path(output_radar, "Function_one_vs_rest_significance.csv"), row.names = FALSE)
write.csv(gsva_cat_sample,  file.path(output_radar, "GSVA_function_score_each_sample.csv"),     row.names = FALSE)
write.csv(gsva_cat_group,   file.path(output_radar, "GSVA_function_score_group_mean.csv"),      row.names = FALSE)
write.csv(baseline_df,      file.path(output_radar, "GSVA_function_baseline.csv"),              row.names = FALSE)

## ############################################################################
## 5. Radar plots (22 axes in manuscript order)
## ############################################################################
stress_order <- c(
  "ROS", "Oxidative phosphorylation", "Mitochondrial function", "Bile acid metabolism",
  "Cholesterol metabolism", "Fatty acid oxidation", "Lipid biosynthesis",
  "Antigen presentation", "Chemokines / Chemotaxis", "Immune cell infiltration",
  "EMT", "TGF-beta signaling", "Fibrosis", "Extracellular matrix remodeling",
  "Apoptosis", "Autophagy", "UPR", "ER stress", "Interferon signaling",
  "NF-kappaB signaling", "TNF signaling", "Inflammatory response")

function_map <- data.frame(Function = stress_order, ID = seq_along(stress_order),
                           stringsAsFactors = FALSE)

dark_colors <- c(
  "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
  "#A65628", "#F781BF", "#999999", "#66C2A5", "#FC8D62", "#8DA0CB",
  "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", "#1B9E77",
  "#D95F02", "#7570B3", "#E7298A", "#66A61E")

color_map <- data.frame(ID = seq_along(stress_order), Function = stress_order,
                        Color = dark_colors[1:length(stress_order)],
                        stringsAsFactors = FALSE)

y_min <- min(c(gsva_cat_group$mean_GSVA, baseline_df$baseline_GSVA), na.rm = TRUE)
y_max <- max(c(gsva_cat_group$mean_GSVA, baseline_df$baseline_GSVA), na.rm = TRUE)
y_pad <- (y_max - y_min) * 0.10; if (is.na(y_pad) || y_pad == 0) y_pad <- 0.1
y_lim <- c(y_min - y_pad, y_max + y_pad)

output_radar2 <- file.path(output, "Radar_2")
if (!dir.exists(output_radar2)) dir.create(output_radar2, recursive = TRUE)

plots_list <- list()
for (g in group_levels) {
  plots_list[[g]] <- make_radar_plot(g, gsva_cat_group, pairwise_res, baseline_df,
                                     y_lim, output_radar2,
                                     function_map, color_map, stress_order)$plot
}

## ############################################################################
## 6. Legend table + standalone legend plot
## ############################################################################
legend_df <- data.frame(ID = seq_along(stress_order), Function = stress_order,
                        Color = dark_colors[1:length(stress_order)],
                        Label = paste0(seq_along(stress_order), ": ", stress_order),
                        stringsAsFactors = FALSE)
write.csv(legend_df, file.path(output_radar2, "Pathway_Legend_Table.csv"), row.names = FALSE)

legend_df_reverse <- legend_df %>% arrange(desc(ID))
legend_plot <- ggplot() +
  geom_point(data = legend_df_reverse, aes(x = 1, y = -ID, color = I(Color)),
             size = 6, shape = 19) +
  geom_text(data = legend_df_reverse, aes(x = 1.3, y = -ID, label = Label),
            hjust = 0, size = 4.5, fontface = "bold", color = "black") +
  xlim(0, 3) + ylim(-23, 0) +
  ggtitle("Pathway Legend (1-22)") +
  theme_void() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16, color = "black"),
        plot.margin = margin(20, 20, 20, 20))

ggsave(file.path(output_radar2, "Pathway_Legend_Simple.png"), legend_plot,
       width = 10, height = 14, dpi = 300)
ggsave(file.path(output_radar2, "Pathway_Legend_Simple.pdf"), legend_plot,
       width = 10, height = 14)

## ############################################################################
## 7. Functional GSVA heatmap
## ############################################################################
draw_function_heatmap(gsva_cat_sample, out_dir = output, group_levels = group_levels,
                      group_colors = group_colors,
                      file_prefix = "Functional_GSVA_ComplexHeatmap_with_KW_star")

## ############################################################################
## 8. Metabolic-reprogramming heatmap (separate GSVA result)
## ############################################################################
load(file.path(ORIGINAL_DIR, "02_GSVA_k4", "02_GSVA_k4_MR_gsva_results.Rdata"))
group_df <- load_groups(file.path(ORIGINAL_DIR, "01_MOVICS", "01_MOVICS_cluster_results.csv"))

mr_dir <- file.path(output, "MR")
if (!dir.exists(mr_dir)) dir.create(mr_dir, recursive = TRUE)

# Reuse the heatmap helper by restructuring MR GSVA into the long category format
gsva_mr_long <- as.data.frame(gsva_results) %>%
  tibble::rownames_to_column("Pathway") %>%
  pivot_longer(cols = -Pathway, names_to = "Sample", values_to = "GSVA") %>%
  left_join(group_df, by = c("Sample" = "sample")) %>%
  mutate(Function = Pathway) %>%
  rename(Function_GSVA = GSVA) %>%
  select(Sample, group, Function, Function_GSVA)

draw_function_heatmap(gsva_mr_long, out_dir = mr_dir, group_levels = group_levels,
                      group_colors = group_colors,
                      file_prefix = "GSVA_Heatmap_with_KW_star",
                      title = "Metabolic Reprogramming GSVA Heatmap by CS Subtype")
