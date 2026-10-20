rm(list = ls()); gc()

## -------- Paths & directories --------
ORIGINAL_DIR <- ""                        # working directory (set as needed)
output       <- file.path(ORIGINAL_DIR, "07_survival")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
setwd(ORIGINAL_DIR)

# Purpose: (1) apply the 40-marker classifier to the GSE15654 cirrhosis cohort,
# test subtype differences in overall survival and HCC-free progression with
# KM + pairwise Cox models, and score Fenton/stress/immune features; (2) use
# XGBoost + Cox to identify which scores drive the poor prognosis of CS3.
# k = 4.

## -------- Libraries --------
library(MOVICS)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(survival)
library(survminer)
library(stringr)
library(MASS)
library(GSVA)
library(GSEABase)
library(IOBR)
library(tibble)

## -------- Settings --------
labels       <- "07_survival_"
palette4     <- c("#3C5488", "#DC0000", "#00A087", "#F39B7F")
palette3     <- c("#3C5488", "#DC0000", "#00A087")

## ############################################################################
## Helper 1: pairwise Cox HR / CI / p across all group pairs
## ############################################################################
pairwise_cox <- function(df_plot, groups) {
  pairs <- combn(groups, 2, simplify = FALSE)
  lapply(pairs, function(pr) {
    df_pair <- df_plot %>% filter(group %in% pr)
    df_pair$group <- factor(df_pair$group, levels = pr)
    cox_uni <- coxph(Surv(time, status) ~ group, data = df_pair)
    summ    <- summary(cox_uni)
    list(pair    = paste(pr[1], "vs", pr[2]),
         hr      = summ$coefficients[, "exp(coef)"],
         ci_low  = summ$conf.int[, "lower .95"],
         ci_high = summ$conf.int[, "upper .95"],
         p_txt   = ifelse(is.na(summ$coefficients[, "Pr(>|z|)"]), "NA",
                          format(signif(summ$coefficients[, "Pr(>|z|)"], 3), scientific = TRUE)))
  })
}

## ############################################################################
## Helper 2: draw a KM plot with pairwise-Cox annotation and save png + pdf
## ############################################################################
plot_km <- function(df_plot, groups, palette, legend_labs, ylab, out_file,
                    annot_x = 15, annot_y = 0.05, width = 6, height = 6) {
  fit <- survfit(Surv(time, status) ~ group, data = df_plot)
  cox_results <- pairwise_cox(df_plot, groups)
  p_note <- paste(sapply(cox_results, function(r) {
    sprintf("%s: HR=%.2f (95%% CI %.2f-%.2f), p=%s",
            r$pair, r$hr, r$ci_low, r$ci_high, r$p_txt)
  }), collapse = "\n")

  p_km <- ggsurvplot(
    fit, data = df_plot,
    risk.table       = TRUE,
    pval             = FALSE,
    conf.int         = FALSE,
    surv.median.line = "none",
    palette          = palette,
    xlab             = "Time (years)",
    ylab             = ylab,
    legend.title     = "Cluster",
    legend.labs      = legend_labs,
    risk.table.height = 0.25,
    ggtheme          = theme_minimal(base_size = 13))

  p_km$plot <- p_km$plot +
    annotate("text", x = annot_x, y = annot_y,
             label = str_wrap(p_note, width = 50),
             parse = FALSE, hjust = 1.02, vjust = 0, size = 3.5, lineheight = 0.9) +
    labs(color = "Cluster", fill = "Cluster") +
    theme(legend.title = element_text(face = "bold"))

  png(paste0(out_file, ".png"), width = width, height = height, res = 300, units = "in")
  print(p_km); dev.off()
  pdf(paste0(out_file, ".pdf"), width = width, height = height)
  print(p_km); dev.off()
  invisible(p_km)
}

## ############################################################################
## Helper 3: read a single-line GMT and return a unique gene vector
## ############################################################################
read_gmt_genes <- function(path, skip_first = FALSE) {
  gmt <- read.csv(path, header = FALSE, sep = "\t")
  if (skip_first) gmt <- gmt[-1, ]
  gene_set <- gmt[, -2]
  rownames(gene_set) <- NULL
  gene_set <- gene_set %>% column_to_rownames("V1") %>% t()
  genes <- unique(as.character(unlist(gene_set)))
  genes[nchar(genes) != 0]
}

## ############################################################################
## Helper 4: XGBoost survival model -> feature importance + risk-stratified KM
## ############################################################################
xgb_survival <- function(factor_df, out_dir, feature_importance_csv = "xgboost_feature_importance.csv") {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  factor_df <- factor_df[complete.cases(factor_df), ]

  X            <- as.matrix(factor_df[, !colnames(factor_df) %in% c("sample", "CS3", "days", "death")])
  feature_names <- colnames(X)
  y_time       <- factor_df$days
  y_status     <- factor_df$death

  set.seed(42)
  dmat <- xgb.DMatrix(X, label = y_time)
  xgb_model <- xgb.train(
    params = list(objective = "survival:cox", eval_metric = "cox-nloglik",
                  max_depth = 4, eta = 0.1, subsample = 0.8, colsample_bytree = 0.8),
    data = dmat, nrounds = 100, verbose = 0, watchlist = list(train = dmat))

  importance_matrix <- xgb.importance(model = xgb_model, feature_names = feature_names) %>%
    arrange(desc(Gain)) %>%
    mutate(Feature = factor(Feature, levels = Feature))

  risk_score <- predict(xgb_model, dmat)
  cox_model  <- coxph(Surv(y_time, y_status) ~ risk_score)
  risk_group <- ifelse(risk_score > median(risk_score), "High", "Low")

  write.csv(importance_matrix, file.path(out_dir, feature_importance_csv), row.names = FALSE)
  save(xgb_model, importance_matrix, risk_score, risk_group, cox_model,
       file = file.path(out_dir, "xgboost_cox_model.RData"))

  p_gain <- ggplot(importance_matrix, aes(x = Gain, y = reorder(Feature, Gain))) +
    geom_bar(stat = "identity", fill = "#2C7FB8", alpha = 0.8) +
    geom_text(aes(label = round(Gain, 4)), hjust = -0.1, size = 3.5) +
    labs(title = "XGBoost Feature Importance (Gain)", x = "Gain", y = "Feature") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10),
          panel.grid.major.y = element_blank())
  ggsave(file.path(out_dir, "xgboost_feature_importance_gain.pdf"), p_gain,
         width = 8, height = max(4, nrow(importance_matrix) * 0.4))
  ggsave(file.path(out_dir, "xgboost_feature_importance_gain.png"), p_gain,
         width = 8, height = max(4, nrow(importance_matrix) * 0.4), dpi = 300)

  fit <- survfit(Surv(y_time, y_status) ~ risk_group)
  p_km <- ggsurvplot(
    fit, data = data.frame(time = y_time, status = y_status, risk_group = risk_group),
    pval = TRUE, pval.method = TRUE, risk.table = TRUE,
    title = "Kaplan-Meier Survival Curve by XGBoost Risk Score",
    xlab = "Days", ylab = "Survival Probability",
    legend.title = "Risk Group", legend.labs = c("High Risk", "Low Risk"),
    palette = c("#E41A1C", "#377EB8"), risk.table.col = "strata",
    ggtheme = theme_minimal())
  ggsave(file.path(out_dir, "KM_curve_by_xgboost_risk_score.pdf"), p_km$plot, width = 8, height = 6)
  ggsave(file.path(out_dir, "KM_curve_by_xgboost_risk_score.png"), p_km$plot, width = 8, height = 6, dpi = 300)

  list(importance = importance_matrix, cox = cox_model, risk_score = risk_score, km = p_km)
}

## ############################################################################
## Helper 5: score gene sets + ESTIMATE + MCPcounter, merge into one table
## ############################################################################
score_features <- function(exp, gene_set, out_dir) {
  gsva_results <- as.data.frame(gsva(gsvaParam(as.matrix(exp), gene_set), verbose = TRUE))
  save(gsva_results, file = file.path(out_dir, paste0(labels, "GSE15654_gsva_results.Rdata")))

  # ESTIMATE / MCPcounter require log2 + de-duplicated symbols
  exp_mat <- as.matrix(exp); storage.mode(exp_mat) <- "numeric"
  exp_mat <- exp_mat[rowSums(is.na(exp_mat)) < ncol(exp_mat), , drop = FALSE]
  exp_mat <- exp_mat[rowSums(exp_mat, na.rm = TRUE) > 0, , drop = FALSE]
  exp_mat <- log2(exp_mat + 1)
  if (any(duplicated(rownames(exp_mat)))) {
    exp_mat <- as.data.frame(exp_mat) %>%
      rownames_to_column("gene") %>%
      group_by(gene) %>%
      summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") %>%
      column_to_rownames("gene") %>% as.matrix()
  }
  rownames(exp_mat) <- trimws(rownames(exp_mat))

  estimate_res <- deconvo_estimate(eset = exp_mat, project = "MASLD", platform = "affymetrix")
  rownames(estimate_res) <- estimate_res$ID
  write.csv(estimate_res, file.path(out_dir, "GSE15654_ESTIMATE_scores.csv"), row.names = FALSE)

  mcp_res <- deconvo_mcpcounter(eset = exp_mat, project = "MASLD")
  rownames(mcp_res) <- estimate_res$ID
  write.csv(mcp_res, file.path(out_dir, "GSE15654_MCPcounter_res.csv"))

  gsva_t   <- as.data.frame(t(gsva_results)); gsva_t$sample   <- rownames(gsva_t)
  est_clean <- estimate_res[, -(1:2), drop = FALSE]; est_clean$sample <- rownames(estimate_res)
  mcp_clean <- mcp_res[, -(1:2), drop = FALSE];      mcp_clean$sample <- rownames(mcp_res)

  common_samples <- Reduce(intersect, list(gsva_t$sample, est_clean$sample, mcp_clean$sample))
  gsva_t    <- gsva_t[gsva_t$sample %in% common_samples, ] %>% arrange(sample)
  est_clean <- est_clean[est_clean$sample %in% common_samples, ] %>% arrange(sample)
  mcp_clean <- mcp_clean[mcp_clean$sample %in% common_samples, ] %>% arrange(sample)

  merged_data <- gsva_t %>%
    left_join(est_clean, by = "sample") %>%
    left_join(mcp_clean, by = "sample") %>%
    column_to_rownames("sample")

  merged_dir <- file.path(out_dir, "merged_data")
  if (!dir.exists(merged_dir)) dir.create(merged_dir, recursive = TRUE)
  write.csv(merged_data, file.path(merged_dir, "merged_GSVA_ESTIMATE_MCP.csv"), row.names = TRUE)
  save(merged_data, file = file.path(merged_dir, "merged_GSVA_ESTIMATE_MCP.RData"))
  merged_data
}

## ############################################################################
## 1. Load GSE15654 cohort + apply the 40-marker classifier
## ############################################################################
marker.up <- readRDS(file.path(ORIGINAL_DIR, "01_MOVICS", "marker.up.rData"))

raw_dir <- file.path(output, "rawdata")
exp   <- read.csv(file.path(raw_dir, "00.rawdata_GSE15654_exp.csv"), row.names = 1)
pdata <- read.csv(file.path(raw_dir, "00.rawdata_GSE15654_pdata.csv"), row.names = 1)
group <- read.csv(file.path(raw_dir, "00.rawdata_GSE15654_group.csv"), row.names = 1)
surv  <- read.csv(file.path(raw_dir, "00.rawdata_GSE15654_sur.csv"), row.names = 1)
surv  <- surv[colnames(exp), , drop = FALSE]

set.seed(123)
GSE.ntp.pred <- runNTP(expr       = exp,
                       templates  = marker.up$templates,
                       scaleFlag  = TRUE,
                       centerFlag = TRUE,
                       doPlot     = TRUE,
                       fig.name   = "NTP HEATMAP FOR GSE15654",
                       width = 8, height = 6)

cluster_p   <- GSE.ntp.pred$ntp.res
cluster_res <- GSE.ntp.pred$clust.res
cluster_res$clust <- paste0("CS", cluster_res$clust)
write.csv(cluster_res, file.path(output, "GSE15654_cluster_res.csv"))
write.csv(cluster_p,   file.path(output, "GSE15654_cluster_p.csv"))

## ############################################################################
## 2. KM: overall survival (all 4 subtypes)
## ############################################################################
df_os <- data.frame(time   = surv$days_to_death / 365,
                    status = surv$death,
                    group  = cluster_res$clust)
plot_km(df_os, groups = c("CS1", "CS2", "CS3", "CS4"), palette = palette4,
        legend_labs = c("CS1", "CS2", "CS3", "CS4"),
        ylab = "Survival probability (OS)",
        out_file = file.path(output, "GSE15654_KM_plot_4groups"))

## ############################################################################
## 3. KM: HCC-free progression (all 4 subtypes)
## ############################################################################
df_hcc <- data.frame(time   = surv$days_to_hcc / 365,
                     status = surv$hcc,
                     group  = cluster_p$prediction)
plot_km(df_hcc, groups = c("CS1", "CS2", "CS3", "CS4"), palette = palette4,
        legend_labs = c("CS1", "CS2", "CS3", "CS4"),
        ylab = "Cumulative survival free of HCC progression",
        out_file = file.path(output, "GSE15654_KM_plot_4groups_Cumulative survival free of HCC progression"))

## ############################################################################
## 4. KM restricted to significant assignments (FDR < 0.05), drop CS4
## ############################################################################
cluster_p_T <- subset(cluster_p, cluster_p$FDR < 0.05)
surv_T      <- surv[rownames(cluster_p_T), , drop = FALSE]

df_os_T <- data.frame(time   = surv_T$days_to_death / 365,
                      status = surv_T$death,
                      group  = cluster_p_T$prediction) %>%
  filter(group != "CS4")
plot_km(df_os_T, groups = c("CS1", "CS2", "CS3"), palette = palette3,
        legend_labs = c("CS1", "CS2", "CS3"),
        ylab = "Survival probability (OS)",
        out_file = file.path(output, "GSE15654_KM_plot_4groups_T"))

df_hcc_T <- data.frame(time   = surv_T$days_to_hcc / 365,
                       status = surv_T$hcc,
                       group  = cluster_p_T$prediction) %>%
  filter(group != "CS4")
plot_km(df_hcc_T, groups = c("CS1", "CS2", "CS3"), palette = palette3,
        legend_labs = c("CS1", "CS2", "CS3"),
        ylab = "Cumulative survival free of HCC progression",
        out_file = file.path(output, "GSE15654_KM_plot_4groups_Cumulative survival free of HCC progression_T"))

## ############################################################################
## 5. Build Fenton-reaction and stress gene sets
## ############################################################################
fenton_root <- file.path(ORIGINAL_DIR, "05_Fenton", "gmt")
stress_root <- file.path(ORIGINAL_DIR, "05_stress", "gmt")
gene_set <- list()

gene_set[["cyFR"]]  <- read_gmt_genes(file.path(fenton_root, "Cytosolic_FR.gmt"))
gene_set[["exFR"]]  <- read_gmt_genes(file.path(fenton_root, "Extracellular_FR.gmt"))
gene_set[["mitFR"]] <- read_gmt_genes(file.path(fenton_root, "Mitochondrial_FR.gmt"))
# Membrane-bound FR: use the first two gene rows of the mitochondrial file
gmt <- read.csv(file.path(fenton_root, "Mitochondrial_FR.gmt"), header = FALSE, sep = "\t")
geneSet <- gmt[, -2]; geneSet <- geneSet[c(1, 2), ]; rownames(geneSet) <- NULL
geneSet <- geneSet %>% column_to_rownames("V1") %>% t()
genes <- unique(as.character(unlist(geneSet))); genes <- genes[nchar(genes) != 0]
gene_set[["mrFR"]] <- genes

gene_set[["OXIDATIVE_STRESS"]] <- read_gmt_genes(
  file.path(stress_root, "oxidative_stress_geneset", "genesets.v2026.1.Hs.gmt"), skip_first = TRUE)
gene_set[["HYPOXIA"]] <- read_gmt_genes(
  file.path(stress_root, "hypoxia_geneset", "genesets.v2026.1.Hs.gmt"))
gene_set[["CHRONIC_INFLAMMATORY"]] <- read_gmt_genes(
  file.path(stress_root, "chronic_inflammation_geneset", "genesets.v2026.1.Hs.gmt"))
gene_set[["ER_STRESS"]] <- read_gmt_genes(
  file.path(stress_root, "er_stress_geneset", "genesets.v2026.1.Hs.gmt"))
gene_set[["INTRACELLULAR_PH_ELEVATION"]] <- read_gmt_genes(
  file.path(stress_root, "stress_geneset", "5-pHi.gmt"))
gene_set[["EXTRACELLULAR_PH_ELEVATION"]] <- read_gmt_genes(
  file.path(stress_root, "stress_geneset", "4-pHe.gmt"))

## ############################################################################
## 6. Score features (GSVA + ESTIMATE + MCPcounter) and merge
## ############################################################################
merged_data <- score_features(exp, gene_set, output)

## ############################################################################
## 7. Build the CS3-driver feature table and run XGBoost + Cox
## ############################################################################
common_samples <- intersect(rownames(cluster_p_T), rownames(merged_data))
merged_data <- merged_data[common_samples, , drop = FALSE]
surv_T      <- surv_T[common_samples, , drop = FALSE]
group_T     <- group[common_samples, , drop = FALSE]

factor_df <- data.frame(
  sample = common_samples,
  CS3    = ifelse(cluster_p_T[common_samples, "prediction"] == "CS3", 1, 0),
  merged_data,
  days   = surv_T$days_to_death,
  death  = surv_T$death,
  stringsAsFactors = FALSE)

# Keep the same reduced feature set used in the original analysis
# (CS3 + 11 feature columns + days + death)
factor_df <- factor_df[, c(1, 2, 3, 4, 5, 7, 8, 9, 10, 11, 13, 25, 26, 27)]
write.csv(factor_df, file.path(output, "factor_df_GSE15654.csv"))

xgb_survival(factor_df, out_dir = file.path(output, "CS3_driver_analysis"))
