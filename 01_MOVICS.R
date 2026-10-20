rm(list = ls()); gc()

## -------- Paths & directories --------
ORIGINAL_DIR <- ""                       # working directory (set as needed)
output       <- file.path(ORIGINAL_DIR, "01_MOVICS")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
setwd(ORIGINAL_DIR)

# Purpose: subtype MASLD samples with MOVICS multi-omics clustering, using
# enzyme-gene expression plus ssGSEA metabolic pathway scores, with k = 4.

library(MOVICS)
library(GSVA)

## -------- Input data --------
# Merged, batch-corrected expression matrix
load(file.path(ORIGINAL_DIR, "00_rawdata_2", "exp_merge_combat_nonneg.rData"))
exp_MASLD <- as.data.frame(res)
range(exp_MASLD)

# Pre-computed ssGSEA metabolic pathway scores
load(file.path(ORIGINAL_DIR, "00_KEGG_get", "ssgsea_01.rData"))

# Enzyme gene list (HumanCyc)
enzyme_list <- read.csv(file.path(ORIGINAL_DIR, "00_rawdata",
                                  "metabolism_subtyping_genes", "gendis_101589_mmc2.csv"))
enzyme_gene <- unique(enzyme_list$Genes)
enzyme_gene <- intersect(enzyme_gene, rownames(exp_MASLD))

exp_enzyme <- exp_MASLD[enzyme_gene, ]
mat        <- as.matrix(exp_enzyme)

# Two-omics input for MOVICS: mRNA expression + metabolic pathway scores
ma.data <- list(mRNA.expr         = as.data.frame(mat),
                metabolism_pathway = as.data.frame(ssgsea_01))

## -------- Cluster-number selection --------
library(factoextra)
library(NbClust)

output <- file.path(ORIGINAL_DIR, "01_MOVICS", "cluster_analysis")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
labels <- "cluster_analysis"

# Within-cluster sum of squares (elbow method)
data <- as.data.frame(mat)
p1 <- fviz_nbclust(data, kmeans, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2, color = "red")
ggsave(file.path(output, paste0(labels, "_elbow_method.pdf")), plot = p1, width = 8, height = 6)
ggsave(file.path(output, paste0(labels, "_elbow_method.png")), plot = p1, width = 8, height = 6, dpi = 300)

# Multi-index cluster-number assessment
nb <- NbClust(data, distance = "euclidean", min.nc = 2, max.nc = 10,
              method = "kmeans", index = "all")

# MOVICS automatic cluster-number selection
set.seed(110)
optk.brca <- getClustNum(data        = ma.data,
                         is.binary   = c(FALSE, FALSE),
                         try.N.clust = 2:8,
                         fig.name    = "CLUSTER NUMBER OF MASLD",
                         fig.path    = output)

## -------- Run nine clustering algorithms --------
# iClusterBayes
set.seed(123)
iClusterBayes.res <- getiClusterBayes(data        = ma.data,
                                      N.clust     = 4,
                                      type        = c("gaussian", "gaussian"),
                                      n.burnin    = 1800,
                                      n.draw      = 1200,
                                      prior.gamma = c(0.5, 0.5),
                                      sdev        = 0.05,
                                      thin        = 3)

# Remaining nine algorithms
set.seed(123)
moic.res.list <- getMOIC(data        = ma.data,
                         methodslist = list("SNF", "CIMLR", "PINSPlus", "NEMO", "COCA",
                                            "MoCluster", "LRAcluster", "ConsensusClustering", "IntNMF"),
                         N.clust     = 4,
                         type        = c("gaussian", "gaussian"))

save(moic.res.list, file = file.path(output, "moic.res.list_01.rda"))

moic.res.list <- append(moic.res.list, list("iClusterBayes" = iClusterBayes.res))
save(moic.res.list, file = file.path(output, "moic.res.list_02.rda"))
load(file.path(output, "moic.res.list_02.rda"))

## -------- Consensus clustering --------
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")
save(cmoic.brca, file = file.path(output, "cmoic.brca.rda"))
load(file.path(output, "cmoic.brca.rda"))

getSilhouette(sil      = cmoic.brca$sil,
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height   = 5.5,
              width    = 5)

## -------- Multi-omics heatmap --------
indata   <- ma.data
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2, 2),
                     centerFlag = c(TRUE, TRUE),
                     scaleFlag  = c(TRUE, TRUE))

mRNA.col   <- c("#00FF00", "#008000", "#000000", "#800000", "#FF0000")
lncRNA.col <- c("#6699CC", "white",   "#FF3C38")
col.list   <- list(mRNA.col, lncRNA.col)

getMoHeatmap(data          = plotdata,
             row.title     = c("mRNA", "metabolism_pathway"),
             is.binary     = c(FALSE, FALSE),
             legend.name   = c("mRNA exp", "metabolism_pathway score"),
             clust.res     = cmoic.brca$clust.res,
             clust.dend    = NULL,
             show.rownames = c(FALSE, FALSE),
             show.colnames = FALSE,
             color         = col.list,
             annCol        = NULL,
             annColors     = NULL,
             width         = 10,
             height        = 5,
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")

## -------- Differential expression & markers: mRNA --------
runDEA(dea.method = "limma",
       expr       = exp_MASLD,
       moic.res   = cmoic.brca,
       prefix     = "MASLD")

marker.up <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "limma",
                       prefix        = "MASLD",
                       dat.path      = getwd(),
                       res.path      = getwd(),
                       p.cutoff      = 0.05,
                       p.adj.cutoff  = 0.05,
                       dirct         = "up",
                       n.marker      = 100,
                       doplot        = TRUE,
                       norm.expr     = exp_MASLD,
                       show_rownames = FALSE,
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")

marker.dn <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "limma",
                       prefix        = "MASLD",
                       dirct         = "down",
                       n.marker      = 50,
                       doplot        = TRUE,
                       norm.expr     = exp_MASLD,
                       fig.name      = "DOWNREGULATED BIOMARKER HEATMAP")

write.csv(marker.up$templates, file.path(output, paste0(labels, "marker_up_exp.csv")))

## -------- Differential expression & markers: metabolic pathways --------
runDEA(dea.method = "limma",
       expr       = as.data.frame(ssgsea_01),
       moic.res   = cmoic.brca,
       prefix     = "MASLD_02")

marker_02.up <- runMarker(moic.res      = cmoic.brca,
                          dea.method    = "limma",
                          prefix        = "MASLD_02",
                          dat.path      = getwd(),
                          res.path      = getwd(),
                          p.cutoff      = 0.05,
                          p.adj.cutoff  = 0.05,
                          dirct         = "up",
                          n.marker      = 5,
                          doplot        = TRUE,
                          norm.expr     = as.data.frame(ssgsea_01),
                          show_rownames = TRUE,
                          fig.name      = "UPREGULATED BIOMARKER HEATMAP02")

marker_02.dn <- runMarker(moic.res      = cmoic.brca,
                          dea.method    = "limma",
                          prefix        = "MASLD_02",
                          dirct         = "down",
                          n.marker      = 5,
                          doplot        = TRUE,
                          norm.expr     = as.data.frame(ssgsea_01),
                          show_rownames = TRUE,
                          fig.name      = "DOWNREGULATED BIOMARKER HEATMAP02")

write.csv(marker_02.up$templates, file.path(output, paste0(labels, "marker_up02_pathway.csv")))

## -------- GSEA enrichment (preliminary) --------
MSIGDB.FILE <- system.file("extdata", "c5.bp.v7.1.symbols.xls", package = "MOVICS", mustWork = TRUE)

gsea.up <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "limma",
                   prefix       = "MASLD",
                   dat.path     = getwd(),
                   res.path     = getwd(),
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = exp_MASLD,
                   dirct        = "up",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "gsva",
                   norm.method  = "mean",
                   fig.name     = "UPREGULATED PATHWAY HEATMAP")

gsea.dn <- runGSEA(moic.res     = cmoic.brca,
                   dea.method   = "limma",
                   prefix       = "MASLD",
                   msigdb.path  = MSIGDB.FILE,
                   norm.expr    = exp_MASLD,
                   dirct        = "down",
                   p.cutoff     = 0.05,
                   p.adj.cutoff = 0.25,
                   gsva.method  = "ssgsea",
                   norm.method  = "median",
                   fig.name     = "DOWNREGULATED PATHWAY HEATMAP")

saveRDS(marker.up, file = file.path(output, "marker.up.rData"))

## -------- Cluster composition vs. clinical groups --------
library(dplyr)
library(ggplot2)
library(scales)
library(patchwork)

group      <- read.csv(file.path(ORIGINAL_DIR, "00_rawdata_2", "00_rawdata_group_merged.csv"), row.names = 1)
group_MASLD <- group[colnames(exp_MASLD), ]
table(group_MASLD$group)
group_MASLD <- group_MASLD %>%
  mutate(group = case_when(
    group %in% c("MASLD", "NAFLD") ~ "MASLD",
    group %in% c("NAFL")  ~ "SS",
    group %in% c("NASH")  ~ "MASH",
    TRUE ~ NA_character_))
table(group_MASLD$group)

cluster_results <- as.data.frame(cmoic.brca$clust.res)
cluster_results <- cluster_results[rownames(group_MASLD), ]
write.csv(cluster_results, file.path(output, paste0(labels, "cluster_results.csv")))

group_MASLD$cluster <- as.character(cluster_results$clust)

plot_dir <- file.path(output, "plot_res")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## -------- Composition plots: SS vs. MASH --------
df <- group_MASLD %>% filter(group %in% c("MASH", "SS"))

# Clusters within each group
comp_by_group <- df %>%
  count(group, cluster) %>%
  group_by(group) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p1 <- ggplot(comp_by_group, aes(x = group, y = prop, fill = factor(cluster))) +
  geom_col() +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Group", y = "Proportion", fill = "Cluster") +
  theme_minimal()
ggsave(file.path(plot_dir, paste0(labels, "proportion_of_clusters_within_each_group.png")), p1, width = 6, height = 6)
ggsave(file.path(plot_dir, paste0(labels, "proportion_of_clusters_within_each_group.pdf")), p1, width = 6, height = 6)

# Groups within each cluster
comp_by_cluster <- df %>%
  count(cluster, group) %>%
  group_by(cluster) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p2 <- ggplot(comp_by_cluster, aes(x = factor(cluster), y = prop, fill = group)) +
  geom_col() +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "Cluster", y = "Proportion", fill = "Group") +
  theme_minimal()
ggsave(file.path(plot_dir, paste0(labels, "proportion_of_groups_within_each_cluster.png")), p2, width = 6, height = 6)
ggsave(file.path(plot_dir, paste0(labels, "proportion_of_groups_within_each_cluster.pdf")), p2, width = 6, height = 6)

# Note: not all cohorts contain both MASH and SS, so clusters do not separate them cleanly.

## -------- Composition plots: NAS / fibrosis --------
pdata       <- read.csv(file.path(ORIGINAL_DIR, "00_rawdata_2", "pdata_merged.csv"), row.names = 1)
pdata_MASLD <- pdata[rownames(group_MASLD), ]
table(pdata_MASLD$age_clean)
table(pdata_MASLD$gender_clean)
table(pdata_MASLD$nas)
table(pdata_MASLD$fibrosis_stage)

pdata_MASLD$cluster <- as.character(cluster_results$clust)
pdata_MASLD$nas     <- as.character(pdata_MASLD$nas)

# NAS grade
nas_df <- subset(pdata_MASLD, nas != "")
comp_by_cluster <- nas_df %>%
  count(nas, cluster) %>%
  group_by(nas) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p1 <- ggplot(comp_by_cluster, aes(x = nas, y = prop, fill = factor(cluster))) +
  geom_col() +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "NAS grade", y = "Proportion", fill = "Cluster") +
  theme_minimal()
ggsave(file.path(plot_dir, paste0(labels, "proportion_of_clusters_within_each_NAS.png")), p1, width = 6, height = 6)
ggsave(file.path(plot_dir, paste0(labels, "proportion_of_clusters_within_each_NAS.pdf")), p1, width = 6, height = 6)

# Fibrosis stage
fibro_df <- subset(pdata_MASLD, fibrosis_stage != "")
comp_by_cluster <- fibro_df %>%
  count(fibrosis_stage, cluster) %>%
  group_by(fibrosis_stage) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p1 <- ggplot(comp_by_cluster, aes(x = fibrosis_stage, y = prop, fill = factor(cluster))) +
  geom_col() +
  geom_text(aes(label = scales::percent(prop, accuracy = 0.1)),
            position = position_stack(vjust = 0.5), size = 3, color = "white") +
  scale_y_continuous(labels = percent_format()) +
  labs(x = "fibrosis_stage", y = "Proportion", fill = "Cluster") +
  theme_minimal()
ggsave(file.path(plot_dir, paste0(labels, "proportion_of_clusters_within_each_fibrosis_stage.png")), p1, width = 6, height = 6)
ggsave(file.path(plot_dir, paste0(labels, "proportion_of_clusters_within_each_fibrosis_stage.pdf")), p1, width = 6, height = 6)

## -------- Marker-number optimization (Kappa scan) --------
output_base <- file.path(output, "marker_optimization")
if (!dir.exists(output_base)) dir.create(output_base, recursive = TRUE)

marker_numbers <- seq(10, 150, by = 10)
n_iterations   <- length(marker_numbers)

kappa_results <- data.frame(
  n_marker     = marker_numbers,
  kappa        = numeric(n_iterations),
  accuracy     = numeric(n_iterations),
  ntp_clusters = numeric(n_iterations),
  stringsAsFactors = FALSE)

for (i in seq_along(marker_numbers)) {
  n_marker <- marker_numbers[i]
  cat("\n", strrep("=", 10), "\n")
  cat(sprintf("Iteration %d/%d: n.marker = %d\n", i, n_iterations, n_marker))
  cat(strrep("=", 10), "\n")

  iter_output <- file.path(output_base, paste0("nmarker_", n_marker))
  if (!dir.exists(iter_output)) dir.create(iter_output, recursive = TRUE)

  tryCatch({
    # Step 1: identify markers
    marker.up <- runMarker(
      moic.res      = cmoic.brca,
      dea.method    = "limma",
      prefix        = "MASLD",
      dat.path      = getwd(),
      res.path      = iter_output,
      p.cutoff      = 0.05,
      p.adj.cutoff  = 0.05,
      dirct         = "up",
      n.marker      = n_marker,
      doplot        = TRUE,
      norm.expr     = exp_MASLD,
      show_rownames = FALSE,
      fig.name      = paste0("UPREGULATED_BIOMARKER_HEATMAP_n", n_marker))

    # Step 2: NTP prediction
    ntp.pred <- runNTP(
      expr       = exp_MASLD,
      templates  = marker.up$templates,
      scaleFlag  = TRUE,
      centerFlag = TRUE,
      doPlot     = TRUE,
      fig.name   = paste0("NTP_HEATMAP_n", n_marker),
      width      = 8,
      height     = 6)

    # Step 3: Kappa consistency
    kappa_result <- runKappa(
      subt1     = cmoic.brca$clust.res$clust,
      subt2     = ntp.pred$clust.res$clust,
      subt1.lab = "CSs",
      subt2.lab = "NTP",
      fig.name  = paste0("CONSISTENCY_HEATMAP_n", n_marker),
      fig.path  = iter_output)

    # Step 4: store
    kappa_results$kappa[i]        <- kappa_result$kappa
    kappa_results$accuracy[i]     <- kappa_result$accuracy
    kappa_results$ntp_clusters[i] <- length(unique(ntp.pred$clust.res$clust))

    cat(sprintf("  Kappa: %.4f\n", kappa_result$kappa))
    cat(sprintf("  Accuracy: %.4f\n", kappa_result$accuracy))
    cat(sprintf("  NTP clusters: %d\n", length(unique(ntp.pred$clust.res$clust))))

    saveRDS(kappa_result, file.path(iter_output, "kappa_result.rds"))
    saveRDS(ntp.pred,     file.path(iter_output, "ntp_pred.rds"))
    saveRDS(marker.up,    file.path(iter_output, "marker_up.rds"))
  }, error = function(e) {
    cat(sprintf("ERROR at n.marker = %d: %s\n", n_marker, e$message))
    kappa_results$kappa[i]        <- NA
    kappa_results$accuracy[i]     <- NA
    kappa_results$ntp_clusters[i] <- NA
  })

  gc()
}

# Kappa scan results (hard-coded from a prior run; 40 markers selected)
plot_data <- na.omit(kappa_results)
plot_data$kappa <- c(0.540, 0.519, 0.542, 0.543, 0.539, 0.530, 0.517,
                     0.507, 0.502, 0.494, 0.493, 0.497, 0.498, 0.486, 0.483)
best_n_marker <- 40

## -------- Final classifier (40 markers) & NTP --------
marker.up <- runMarker(moic.res      = cmoic.brca,
                       dea.method    = "limma",
                       prefix        = "MASLD",
                       dat.path      = getwd(),
                       res.path      = getwd(),
                       p.cutoff      = 0.05,
                       p.adj.cutoff  = 0.05,
                       dirct         = "up",
                       n.marker      = 40,
                       doplot        = TRUE,
                       norm.expr     = exp_MASLD,
                       show_rownames = FALSE,
                       fig.name      = "UPREGULATED BIOMARKER HEATMAP")
saveRDS(marker.up, file = file.path(output, "marker.up.rData"))

ntp.final <- runNTP(
  expr       = exp_MASLD,
  templates  = marker.up$templates,
  scaleFlag  = TRUE,
  centerFlag = TRUE,
  doPlot     = TRUE,
  fig.name   = paste0("NTP_HEATMAP_OPTIMAL_n", best_n_marker),
  width      = 8,
  height     = 6)
