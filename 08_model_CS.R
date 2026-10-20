rm(list = ls()); gc()

## -------- Paths & directories --------
ORIGINAL_DIR <- ""                        # working directory (set as needed)
output       <- file.path(ORIGINAL_DIR, "08_model_CS")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
setwd(ORIGINAL_DIR)

# Purpose: build fast classifiers for the CS3 subtype (and, for comparison, all
# four subtypes) from unsupervised-clustering markers; validate AUC, calibration,
# and SHAP-based feature importance in held-out data. Also characterize ADH-family
# expression and the ADH signature score across subtypes and validation cohorts.

## -------- Libraries --------
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(patchwork)
library(xgboost)
library(SHAPforxgboost)
library(pROC)
library(ggpubr)
library(rstatix)
library(grid)
library(GSVA)
library(msigdbr)
library(fgsea)
library(ggridges)

## -------- Settings --------
labels       <- "08_model_CS_"
cols         <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F")
adh_genes    <- c("ADH1A", "ADH1B", "ADH1C", "ADH4", "ADH6")
gse_list     <- c("GSE193066", "GSE213621", "GSE83452", "GSE164760")

## -------- Load expression + subtype labels --------
load(file.path(ORIGINAL_DIR, "00_rawdata_2", "exp_merge_combat_nonneg.rData"))
exp_MASLD <- as.data.frame(res)

group <- read.csv(file.path(ORIGINAL_DIR, "01_MOVICS", "01_MOVICS_cluster_results.csv"), row.names = 1)
group$clust <- paste0("CS", group$clust)
table(group$clust)

## ============================================================================
## Helper 1: build a binary XGBoost classifier for one target subtype
## ============================================================================
build_subtype_classifier <- function(exp_df, group_df, target_subtype,
                                     seed = 123, nrounds = 300) {
  common_samples <- intersect(colnames(exp_df), rownames(group_df))
  exp_df   <- exp_df[, common_samples, drop = FALSE]
  group_sub <- group_df[common_samples, , drop = FALSE]
  y <- ifelse(group_sub$clust == target_subtype, 1, 0)

  x <- t(as.matrix(exp_df)); storage.mode(x) <- "numeric"
  x <- x[complete.cases(x), , drop = FALSE]
  y <- y[complete.cases(x)]
  x <- x[, apply(x, 2, sd, na.rm = TRUE) > 0, drop = FALSE]
  keep_idx <- complete.cases(x) & !is.na(y)
  x <- x[keep_idx, , drop = FALSE]; y <- y[keep_idx]
  stopifnot(nrow(x) == length(y))

  set.seed(seed)
  train_idx <- sample(seq_len(nrow(x)), size = floor(0.7 * nrow(x)))
  x_train <- x[train_idx, , drop = FALSE]; y_train <- y[train_idx]
  x_test  <- x[-train_idx, , drop = FALSE]; y_test <- y[-train_idx]

  params <- list(objective = "binary:logistic", eval_metric = "auc",
                 booster = "gbtree", eta = 0.05, max_depth = 3,
                 min_child_weight = 1, subsample = 0.8, colsample_bytree = 0.8,
                 gamma = 0, alpha = 0.1, lambda = 1)
  params$scale_pos_weight <- sum(y_train == 0) / sum(y_train == 1)

  set.seed(seed)
  model <- xgb.train(params = params, data = xgb.DMatrix(x_train, label = y_train),
                     nrounds = nrounds, verbose = 0)

  test_auc <- as.numeric(pROC::auc(pROC::roc(y_test, predict(model, xgb.DMatrix(x_test)), quiet = TRUE)))

  imp <- xgb.importance(feature_names = colnames(x_train), model = model)
  imp <- imp[order(-imp$Gain), , drop = FALSE]
  selected_genes <- unique(imp$Feature[imp$Gain > 0.003])
  if (length(selected_genes) < 5) selected_genes <- head(imp$Feature, 20)

  x_train_sel <- x_train[, selected_genes, drop = FALSE]
  x_test_sel  <- x_test[,  selected_genes, drop = FALSE]

  set.seed(seed)
  final_model <- xgb.train(params = params,
                           data = xgb.DMatrix(x_train_sel, label = y_train),
                           nrounds = nrounds, verbose = 0)
  final_pred <- predict(final_model, xgb.DMatrix(x_test_sel))
  final_auc  <- as.numeric(pROC::auc(pROC::roc(y_test, final_pred, quiet = TRUE)))

  shap_res <- shap.values(xgb_model = final_model, X_train = x_train_sel)
  shap_long <- shap.prep(shap_contrib = shap_res$shap_score, X_train = x_train_sel)
  shap_summary <- data.frame(gene = names(shap_res$mean_shap_score),
                             mean_shap_score = as.numeric(shap_res$mean_shap_score),
                             row.names = NULL, stringsAsFactors = FALSE) %>%
    arrange(desc(mean_shap_score))

  list(
    model = final_model, selected_genes = selected_genes,
    shap_long = shap_long, shap_summary = shap_summary,
    x_train_sel = x_train_sel, x_test_sel = x_test_sel,
    y_train = y_train, y_test = y_test,
    test_auc_all = test_auc, test_auc_sel = final_auc)
}

## ============================================================================
## Helper 2: calibration plot (training + testing)
## ============================================================================
plot_calibration <- function(pred_train, y_train, pred_test, y_test) {
  r2_train <- cor(pred_train, y_train)^2
  r2_test  <- cor(pred_test,  y_test)^2
  train_df <- data.frame(set = "Training", pred = pred_train, obs = y_train)
  test_df  <- data.frame(set = "Testing",  pred = pred_test,  obs = y_test)

  plot_one <- function(df, title_text, r2) {
    ggplot(df, aes(x = pred, y = obs)) +
      geom_point(alpha = 0.25, size = 1.2, color = "grey50") +
      geom_smooth(method = "loess", se = FALSE, color = "#D7301F", linewidth = 1) +
      geom_abline(intercept = 0, slope = 1, linetype = 2, color = "grey40") +
      coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
      labs(x = "Predicted probability", y = "Observed outcome",
           title = title_text, subtitle = paste0("R² = ", round(r2, 3))) +
      theme_bw(base_size = 14) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            plot.subtitle = element_text(hjust = 0.5))
  }
  list(train = plot_one(train_df, "Training set calibration", r2_train),
       test  = plot_one(test_df,  "Testing set calibration",  r2_test),
       r2_train = r2_train, r2_test = r2_test)
}

## ============================================================================
## Helper 3: per-gene boxplot by subtype with pairwise Wilcoxon brackets
## ============================================================================
plot_gene_box <- function(df, gene, title = gene, ylab = NULL,
                          cols = c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F"),
                          out_dir, prefix = "", ylim_mult = 2.4) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  group_levels <- unique(df$group)
  if (length(cols) < length(group_levels)) cols <- rep(cols, length.out = length(group_levels))
  names(cols) <- group_levels

  fml  <- as.formula(paste(gene, "~ group"))
  kw_p <- kruskal.test(fml, data = df)$p.value
  pairwise_res <- df %>%
    pairwise_wilcox_test(formula = fml, p.adjust.method = "BH") %>%
    filter(p.adj < 0.05) %>%
    add_significance("p.adj") %>%
    add_xy_position(x = "group", step.increase = 0.18)

  p <- ggplot(df, aes(x = group, y = .data[[gene]], fill = group, colour = group)) +
    geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.5,
                 colour = "black", linewidth = 0.4) +
    geom_jitter(position = position_jitter(width = 0.15, height = 0), size = 1.8, alpha = 0.9) +
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    theme_minimal(base_size = 15) +
    labs(x = "", y = ylab %||% paste0("Relative ", gene, " Expression"),
         title = title, subtitle = paste0("Kruskal-Wallis p = ", signif(kw_p, 3))) +
    theme(panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.8),
          axis.title.x    = element_text(size = 16, face = "bold"),
          axis.title.y    = element_text(size = 14),
          axis.text.x     = element_text(size = 14, face = "bold"),
          axis.text.y     = element_text(size = 14),
          legend.position = "none",
          plot.title      = element_text(hjust = 0.5, size = 16, face = "bold"),
          plot.subtitle   = element_text(hjust = 0, size = 12)) +
    coord_cartesian(clip = "off",
                    ylim = c(min(df[[gene]], na.rm = TRUE),
                             max(df[[gene]], na.rm = TRUE) * ylim_mult))

  if (nrow(pairwise_res) > 0) {
    p <- p + stat_pvalue_manual(pairwise_res, label = "p.adj.signif",
                                tip.length = 0.01, hide.ns = TRUE, size = 5)
  }
  ggsave(file.path(out_dir, paste0(prefix, gene, "_plot.png")), p, width = 4, height = 4, dpi = 300)
  ggsave(file.path(out_dir, paste0(prefix, gene, "_plot.pdf")), p, width = 4, height = 4)
  invisible(p)
}

## ============================================================================
## Helper 4: run a full subtype classifier and persist all outputs
## ============================================================================
run_and_save_subtype <- function(exp_df, group_df, target_subtype, out_dir, seed = 123) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  res <- build_subtype_classifier(exp_df, group_df, target_subtype, seed = seed)

  cal <- plot_calibration(
    pred_train = predict(res$model, xgb.DMatrix(res$x_train_sel)),
    y_train    = res$y_train,
    pred_test  = predict(res$model, xgb.DMatrix(res$x_test_sel)),
    y_test     = res$y_test)

  write.csv(data.frame(gene = res$selected_genes),
            file.path(out_dir, paste0(target_subtype, "_selected_genes.csv")), row.names = FALSE)
  write.csv(res$shap_summary,
            file.path(out_dir, paste0(target_subtype, "_SHAP_summary.csv")), row.names = FALSE)
  write.csv(data.frame(subtype = target_subtype,
                       test_auc_all_features = res$test_auc_all,
                       test_auc_selected_genes = res$test_auc_sel,
                       selected_genes_count = length(res$selected_genes),
                       r2_train = cal$r2_train, r2_test = cal$r2_test),
            file.path(out_dir, paste0(target_subtype, "_model_performance.csv")), row.names = FALSE)

  pdf(file.path(out_dir, paste0(target_subtype, "_SHAP_summary.pdf")), width = 8, height = 10)
  print(shap.plot.summary(res$shap_long)); dev.off()
  png(file.path(out_dir, paste0(target_subtype, "_SHAP_summary.png")),
      width = 8, height = 10, units = "in", res = 300)
  print(shap.plot.summary(res$shap_long)); dev.off()

  p_cal <- cal$train + cal$test
  ggsave(file.path(out_dir, paste0(target_subtype, "_Calibration_curve_train_test.png")),
         p_cal, width = 8, height = 4.5, dpi = 300)
  ggsave(file.path(out_dir, paste0(target_subtype, "_Calibration_curve_train_test.pdf")),
         p_cal, width = 8, height = 4.5)
  res
}

## ############################################################################
## 1. Volcano plot of CS3-vs-others DEGs (limma)
## ############################################################################
limma_results <- read.csv(
  file.path(ORIGINAL_DIR, "01_MOVICS", "consensusMOIC_MASLD_limma_test_result.CS3_vs_Others.txt"),
  sep = "\t", row.names = 1)
nrDEG <- limma_results
nrDEG$change <- ifelse(nrDEG$padj < 0.05 & abs(nrDEG$log2fc) >= 0.5,
                       ifelse(nrDEG$log2fc > 0, "Up", "Down"), "NS")
write.csv(nrDEG, file.path(output, "limma_results.csv"))

p1 <- ggplot(nrDEG, aes(x = log2fc, y = -log10(padj), colour = change)) +
  geom_point(alpha = 0.4, size = 2) +
  scale_color_manual(values = c("#4DBBD5", "#d2dae2", "#E64B35")) +
  geom_vline(xintercept = 0,   lty = 4, col = "black",  lwd = 0.8) +
  geom_vline(xintercept = 0.5, lty = 2, col = "grey50", lwd = 0.4) +
  geom_vline(xintercept = -0.5,lty = 2, col = "grey50", lwd = 0.4) +
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +
  labs(x = expression(log[2]~FC), y = expression(-log[10]~italic(p))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "right", legend.text = element_text(size = 12),
        legend.title = element_blank())
ggsave(file.path(output, "DEG_valcano.png"), p1, width = 5, height = 5, units = "in")
ggsave(file.path(output, "DEG_valcano.pdf"), p1, width = 5, height = 5, units = "in")

## ############################################################################
## 2. CS3-only classifier on all genes
## ############################################################################
run_and_save_subtype(exp_MASLD, group, "CS3",
                     out_dir = file.path(output, "CS3_xgboost_SHAP"), seed = 123)

## ############################################################################
## 3. All-subtype classifiers (CS1-CS4) with summary table
## ############################################################################
for (s in c("CS1", "CS2", "CS3", "CS4")) {
  message("Processing: ", s)
  run_and_save_subtype(exp_MASLD, group, s,
                       out_dir = file.path(output, paste0(s, "_xgboost_SHAP")), seed = 67)
}

perf_files <- list.files(output, pattern = "_model_performance.csv", recursive = TRUE, full.names = TRUE)
all_perf <- do.call(rbind, lapply(perf_files, read.csv))
write.csv(all_perf, file.path(output, "all_subtypes_model_performance_summary.csv"), row.names = FALSE)

## ############################################################################
## 4. ADH-family expression across subtypes (discovery)
## ############################################################################
ADH_data <- as.data.frame(t(exp_MASLD[adh_genes, ]))
ADH_data <- as.data.frame(scale(ADH_data))
ADH_data$group <- group$clust

adh_plot_dir <- file.path(output, "ADH_plot")
for (g in adh_genes) {
  plot_gene_box(ADH_data, g, cols = cols, out_dir = adh_plot_dir)
}

## ############################################################################
## 5. ADH-family expression across subtypes (validation cohorts)
## ############################################################################
for (GSE in gse_list) {
  exp  <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(GSE, "_exp.csv")), row.names = 1)
  grp  <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(GSE, "_group.csv")), row.names = 1)
  cg   <- read.csv(file.path(ORIGINAL_DIR, "04_test", GSE, paste0(labels, GSE, "cluster_res.csv")), row.names = 1)
  cg$clust <- paste0("CS", cg$clust)
  cg <- cg[rownames(grp), ]; grp$clust <- cg$clust

  ADH_data <- na.omit(exp[adh_genes, ])
  ADH_data <- as.data.frame(t(ADH_data))
  ADH_data <- as.data.frame(scale(ADH_data))
  ADH_data$group <- grp$clust

  out_dir <- file.path(output, "ADH_plot", GSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  for (g in adh_genes) plot_gene_box(ADH_data, g, cols = cols, out_dir = out_dir, prefix = GSE)
}

## ############################################################################
## 6. ADH signature score (GSVA) + boxplot across subtypes
## ############################################################################
exp_mat <- as.matrix(exp_MASLD); storage.mode(exp_mat) <- "numeric"
if (any(duplicated(rownames(exp_mat)))) {
  exp_df <- as.data.frame(exp_mat); exp_df$gene <- rownames(exp_mat)
  exp_df <- aggregate(. ~ gene, data = exp_df, FUN = mean, na.rm = TRUE)
  rownames(exp_df) <- exp_df$gene; exp_df$gene <- NULL
  exp_mat <- as.matrix(exp_df)
}
adh_set <- list(ADH_signature = intersect(adh_genes, rownames(exp_mat)))
if (length(adh_set[[1]]) < 1) stop("No ADH genes matched the expression matrix rownames.")

adh_gsva  <- gsva(gsvaParam(exp_mat, adh_set), verbose = TRUE)
adh_score <- as.data.frame(t(adh_gsva))
colnames(adh_score) <- "ADH_signature_score"
adh_score$Sample <- rownames(adh_score)
adh_score$group  <- group[rownames(adh_score), "clust"]

adh_gsva_dir <- file.path(output, "ADH_GSVA_plot")
plot_gene_box(adh_score, "ADH_signature_score", title = "",
              ylab = "ADH_signature_score", cols = cols,
              out_dir = adh_gsva_dir, ylim_mult = 1.8)
write.csv(adh_score, file.path(adh_gsva_dir, "adh_gsva_results.csv"))

## ############################################################################
## 7. Pre-ranked GSEA of the ADH signature score (C5 BP)
## ############################################################################
target_vec <- as.data.frame(t(adh_gsva)); mode(target_vec) <- "numeric"
cors <- apply(exp_MASLD, 1, function(x) {
  suppressWarnings(cor(x, target_vec, method = "spearman", use = "pairwise.complete.obs"))
})
ranked <- sort(cors, decreasing = TRUE)

pathways_list <- list()
c5bp <- tryCatch(msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP"),
                 error = function(e) NULL)
save(c5bp, file = file.path(output, "c5bp.Rdata"))
if (!is.null(c5bp)) {
  c5_sets <- split(c5bp$gene_symbol, c5bp$gs_name)
  names(c5_sets) <- paste0("C5_BP|", names(c5_sets))
  pathways_list <- c(pathways_list, c5_sets)
}

set.seed(123)
fgsea_res <- fgsea(pathways = pathways_list, stats = ranked,
                   minSize = 15, maxSize = 500, nperm = 10000) %>%
  as.data.frame() %>% arrange(padj, pval, desc(NES))

safe_fgsea <- fgsea_res
for (col in names(safe_fgsea)[sapply(safe_fgsea, is.list)]) {
  safe_fgsea[[col]] <- sapply(safe_fgsea[[col]], function(x) {
    if (is.null(x) || length(x) == 0) NA_character_ else paste(as.character(x), collapse = ";")
  }, USE.NAMES = FALSE)
}
write.csv(as.data.frame(safe_fgsea, stringsAsFactors = FALSE),
          file.path(adh_gsva_dir, paste0(labels, "adh_set_GSEA_results01.csv")), row.names = FALSE)

sig <- fgsea_res %>% filter(pval < 0.01)
if (nrow(sig) == 0) stop("No pathways with pval < 0.01")

selected <- bind_rows(sig %>% arrange(desc(NES)) %>% head(5),
                      sig %>% arrange(NES) %>% head(5)) %>%
  distinct(pathway, .keep_all = TRUE)
selected_paths  <- selected$pathway
display_names   <- sub("^C5_BP\\|GOBP_", "", selected_paths)
display_names   <- sub("^C5_BP\\|", "", display_names)
display_names   <- str_replace_all(str_to_lower(display_names), "_", " ")
names(display_names) <- selected_paths

eps <- 1e-300
neglog10pval <- -log10(selected$pval + eps); names(neglog10pval) <- selected$pathway

calc_running_es <- function(stats, geneset) {
  N <- length(stats); selected_idx <- which(names(stats) %in% geneset); Nh <- length(selected_idx)
  if (Nh == 0) return(rep(0, N))
  weights <- abs(stats)^1; hit <- integer(N); hit[selected_idx] <- 1
  Phit  <- cumsum(hit * weights) / sum(weights[hit == 1])
  Pmiss <- cumsum((1 - hit) / (N - Nh))
  Phit - Pmiss
}

df <- bind_rows(lapply(selected_paths, function(pw) {
  es <- calc_running_es(ranked, pathways_list[[pw]])
  data.frame(rank = seq_along(ranked), ES = es, ESpos = es - min(es),
             pathway = pw, pathway_display = display_names[pw],
             NES = selected$NES[match(pw, selected$pathway)],
             neglog10pval = neglog10pval[pw], stringsAsFactors = FALSE)
}))
pos_order <- selected %>% arrange(desc(NES)) %>% pull(pathway)
df$pathway_display <- factor(df$pathway_display, levels = display_names[pos_order])
df$pathway_display_wrapped <- str_wrap(df$pathway_display, width = 40)

p <- ggplot(df, aes(x = rank, y = pathway_display_wrapped, height = ESpos, fill = neglog10pval)) +
  geom_ridgeline(stat = "identity", scale = 1, colour = "black", linewidth = 0.2) +
  scale_fill_gradient(low = "#4DBBD5", high = "#E64B35",
                      name = expression(-log[10](italic(p) - value))) +
  labs(x = "Gene rank (by correlation)", y = "",
       title = bquote(italic("ADH_signature_score")~"GSEA")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10, face = "bold"),
        plot.title  = element_text(hjust = 0.5),
        strip.text  = element_text(size = 15, face = "bold"))
ggsave(file.path(adh_gsva_dir, paste0(labels, "adh_set_GSEA_results01.png")), p, width = 6.5, height = 5)
ggsave(file.path(adh_gsva_dir, paste0(labels, "adh_set_GSEA_results01.pdf")), p, width = 6.5, height = 5)
