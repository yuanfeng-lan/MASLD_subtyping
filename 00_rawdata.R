rm(list = ls()); gc()

## -------- Paths & directories --------
ORIGINAL_DIR <- ""                      # working directory (set as needed)
output       <- file.path(ORIGINAL_DIR, "00_rawdata_2")
if (!dir.exists(output)) dir.create(output, recursive = TRUE)
setwd(ORIGINAL_DIR)

# Purpose: integrate high-throughput sequencing samples and normalize them jointly.
file_path <- file.path("", "results", "MASLD")   # raw per-dataset output root
out_dir   <- file.path(output, "MASLD")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

ann   <- read.csv(file.path(output, "gene_annotation_v36.csv"))
geo   <- c("GSE126848", "GSE130970", "GSE135251", "GSE162694", "GSE167523",
           "GSE174478", "GSE192959", "GSE193066", "GSE213621", "GSE225740",
           "GSE240729", "GSE249997", "GSE267031", "GSE276114")

library(edgeR)
library(dplyr)

# Helper: CPM normalization without extra scaling
cpm_norm <- function(counts) {
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "none")
  cpm(dge, normalized.lib.sizes = TRUE)
}

# Helper: write the three standardized output tables for one dataset
save_dataset <- function(GSE, exp, group, pdata) {
  write.csv(exp,   file.path(out_dir, paste0(GSE, "_exp.csv")))
  write.csv(group, file.path(out_dir, paste0(GSE, "_group.csv")))
  write.csv(pdata, file.path(out_dir, paste0(GSE, "_pdata.csv")))
}

## -------- GSE126848 --------
GSE <- "GSE126848"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
colnames(exp) <- rownames(pdata)

cpm_data <- cpm_norm(exp)
table(group$group)
group <- group %>%
  mutate(group = case_when(
    group %in% c("healthy") ~ "Healthy",
    group %in% c("obese")   ~ "Obese",
    group %in% c("NAFLD")   ~ "NAFL",
    group %in% c("NASH")    ~ "NASH",
    TRUE ~ NA_character_))
group <- subset(group, group %in% c("NAFL", "NASH"))
table(group$group)
exp   <- exp[, rownames(group)]
pdata <- pdata[rownames(group), ]
save_dataset(GSE, cpm_data, group, pdata)

## -------- GSE130970 --------
GSE <- "GSE130970"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
colnames(exp) <- rownames(pdata)

table(group$group)
group$group <- "MASLD"
save_dataset(GSE, exp, group, pdata)

## -------- GSE135251 --------
GSE <- "GSE135251"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_counts_data.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
cpm_data <- cpm_norm(exp)

group$group <- pdata$characteristics_ch1.2
table(group$group)
group <- group %>%
  mutate(group = case_when(
    group %in% c("group in paper: control") ~ "Healthy",
    group %in% c("undefined") ~ "Obese",
    group %in% c("group in paper: NAFL") ~ "NAFL",
    group %in% c("group in paper: NASH_F0-F1", "group in paper: NASH_F2",
                 "group in paper: NASH_F3", "group in paper: NASH_F4") ~ "NASH",
    TRUE ~ NA_character_))
table(group$group)
group <- subset(group, group %in% c("NAFL", "NASH"))
table(group$group)
exp   <- exp[, rownames(group)]
pdata <- pdata[rownames(group), ]
save_dataset(GSE, cpm_data, group, pdata)

## -------- GSE162694 --------
GSE <- "GSE162694"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
colnames(exp) <- rownames(pdata)

cpm_data <- cpm_norm(exp)
group$group <- pdata$characteristics_ch1.3
table(group$group)
group$group <- ifelse(group$group == "fibrosis stage: normal liver histology", "Healthy", "NASH")
table(group$group)
group <- subset(group, group == "NASH")
table(group$group)
save_dataset(GSE, cpm_data, group, pdata)

## -------- GSE167523 --------
GSE <- "GSE167523"
exp    <- read.table(file.path(file_path, GSE, paste0(GSE, "_Raw_gene_counts_matrix.txt")),
                     sep = "\t", header = TRUE)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
exp_filter <- aggregate(exp[2:ncol(exp)], list(exp$Gene), FUN = "mean")
rownames(exp_filter) <- exp_filter$Group.1
exp_filter <- exp_filter[, -1]
colnames(exp_filter) <- rownames(pdata)

cpm_data <- cpm_norm(exp_filter)
group <- data.frame(row.names = rownames(pdata), group = pdata$disease.subtype.ch1)
table(group$group)
save_dataset(GSE, cpm_data, group, pdata)

## -------- GSE174478 --------
GSE <- "GSE174478"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)

sample <- pdata[, c(1, 2)]
sample$title <- sub("Liver transcriptome of ", "", sample$title)
intersect(sample$title, colnames(exp))
annotation_map <- setNames(sample$geo_accession, sample$title)
colnames(exp) <- ifelse(colnames(exp) %in% names(annotation_map),
                        annotation_map[colnames(exp)], colnames(exp))
exp <- exp[, rownames(pdata)]
group$group <- "NAFLD"
table(group$group)

cpm_data <- cpm_norm(exp)
save_dataset(GSE, cpm_data, group, pdata)

## -------- GSE192959 --------
GSE <- "GSE192959"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
exp <- exp[, -1]
annotation_map <- setNames(pdata$geo_accession, pdata$title)
colnames(exp) <- ifelse(colnames(exp) %in% names(annotation_map),
                        annotation_map[colnames(exp)], colnames(exp))
group$group <- "NAFLD"
save_dataset(GSE, exp, group, pdata)

## -------- GSE193066 --------
GSE <- "GSE193066"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp_1.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
exp <- exp[, -1]
annotation_map <- setNames(pdata$geo_accession, pdata$title)
colnames(exp) <- ifelse(colnames(exp) %in% names(annotation_map),
                        annotation_map[colnames(exp)], colnames(exp))
total_samples <- ncol(exp)
group$group <- "NAFLD"
save_dataset(GSE, exp, group, pdata)

## -------- GSE213621 --------
GSE <- "GSE213621"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)

sample <- pdata[, 1:2]
sample$title <- sub("human liver ", "", sample$title)
sample$title <- sub(" ", "_", sample$title)
sample$title <- sub("s", "S", sample$title)
annotation_map <- setNames(sample$geo_accession, sample$title)
colnames(exp) <- ifelse(colnames(exp) %in% names(annotation_map),
                        annotation_map[colnames(exp)], colnames(exp))
# Manual fix for a single mislabeled column (column index shifted after mapping)
new_colnames <- c("GSM6589941", colnames(exp))
new_colnames <- new_colnames[-369]
colnames(exp) <- new_colnames
exp <- exp[, rownames(pdata)]
table(group$group)
group <- group %>%
  mutate(group = case_when(
    group %in% c("Control") ~ "Healthy",
    group %in% c("undefined") ~ "Obese",
    group %in% c("F0F1", "F2", "F3F4") ~ "NAFLD",
    group %in% c("NASH liver") ~ "NASH",
    group %in% c("NASH-HCC tumor") ~ "HCC",
    TRUE ~ NA_character_))
table(group$group)
group <- subset(group, group == "NAFLD")
exp   <- exp[, rownames(group)]
pdata <- pdata[rownames(group), ]
save_dataset(GSE, exp, group, pdata)

## -------- GSE225740 --------
GSE <- "GSE225740"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
colnames(exp) <- rownames(pdata)
cpm_data <- cpm_norm(exp)
group$group <- "NAFLD"
table(group$group)
save_dataset(GSE, cpm_data, group, pdata)

## -------- GSE240729 --------
GSE <- "GSE240729"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp_1.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
sample <- pdata[, 1:2]
sample$title <- paste("s", sample$title)
sample$title <- sub(" ", "", sample$title)
sample$title <- gsub("-", ".", sample$title)
annotation_map <- setNames(sample$geo_accession, sample$title)
colnames(exp) <- ifelse(colnames(exp) %in% names(annotation_map),
                        annotation_map[colnames(exp)], colnames(exp))
cpm_data <- cpm_norm(exp)
group$group <- "NAFLD"
save_dataset(GSE, cpm_data, group, pdata)

## -------- GSE249997 --------
GSE <- "GSE249997"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp_1.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
library(stringr)
sample <- pdata[, 1:2]
sample$title <- str_extract(sample$title, "(?<=\\[).*?(?=\\])")
sample$title <- gsub("-", ".", sample$title)
sample$title <- paste("X", sample$title)
sample$title <- sub(" ", "", sample$title)
annotation_map <- setNames(sample$geo_accession, sample$title)
colnames(exp) <- ifelse(colnames(exp) %in% names(annotation_map),
                        annotation_map[colnames(exp)], colnames(exp))
cpm_data <- cpm_norm(exp)
group$group <- "NAFLD"
save_dataset(GSE, cpm_data, group, pdata)

## -------- GSE267031 --------
GSE <- "GSE267031"
exp    <- read.csv(file.path(file_path, GSE, paste0(GSE, "_exp_1.csv")), row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
colnames(exp) <- rownames(pdata)
table(group$group)
group$group <- "NAFLD"
save_dataset(GSE, exp, group, pdata)

## -------- GSE276114 --------
GSE <- "GSE276114"
exp    <- read.table(file.path(file_path, GSE, paste0(GSE, "_raw.count.txt")),
                     sep = "\t", header = TRUE, row.names = 1)
pdata  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_pdata.csv")), row.names = 1)
group  <- read.csv(file.path(file_path, GSE, paste0(GSE, "_group.csv")), row.names = 1)
colnames(exp) <- rownames(pdata)
cpm_data <- cpm_norm(exp)
table(pdata$characteristics_ch1.2)
pdata <- subset(pdata, characteristics_ch1.2 == "disease: MASLD")
exp   <- exp[, rownames(pdata)]
group <- group[rownames(pdata), , drop = FALSE]
group$group <- "NAFLD"
save_dataset(GSE, cpm_data, group, pdata)
