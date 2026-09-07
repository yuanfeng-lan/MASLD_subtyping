#' MASLD Subtype Classification Using 40-Gene Signature
#'
#' This function assigns MASLD samples to one of four metabolic subtypes
#' (CS1-CS4) based on a 40-gene signature trained on 955 samples across
#' 11 GEO cohorts.
#'
#' @param expr_matrix A normalized gene expression matrix with genes as rows
#'        and samples as columns. Row names must be HGNC gene symbols.
#' @param method Classification method. One of "NTP" (Nearest Template
#'        Prediction, default) or "PAM" (Prediction Analysis for Microarrays).
#' @param output_dir Optional directory path for saving output files.
#'        If NULL, no files are written.
#' @param do_plot Logical. Whether to generate a heatmap of classification
#'        results. Default is TRUE.
#'
#' @return A list containing:
#'   \item{cluster_results}{Data frame with sample IDs and assigned subtypes}
#'   \item{heatmap}{Heatmap plot object (if do_plot = TRUE)}
#'
#' @examples
#' \dontrun{
#' # Load expression data
#' expr <- read.csv("path/to/expression_matrix.csv", row.names = 1)
#' 
#' # Classify samples
#' results <- classify_masld_subtype(expr_matrix = expr)
#' 
#' # View subtype distribution
#' table(results$cluster_results$clust)
#' }
#'
#' @export
classify_masld_subtype <- function(expr_matrix,
                                   method = c("NTP", "PAM"),
                                   output_dir = NULL,
                                   do_plot = TRUE) {
  
  # ---- 1. Input Validation ----
  method <- match.arg(method)
  
  if (!is.matrix(expr_matrix) && !is.data.frame(expr_matrix)) {
    stop("expr_matrix must be a matrix or data.frame")
  }
  
  if (is.null(rownames(expr_matrix))) {
    stop("expr_matrix must have gene symbols as row names")
  }
  
  message("Input expression matrix: ", nrow(expr_matrix), " genes, ",
          ncol(expr_matrix), " samples")
  
  # ---- 2. Load Required Resources ----
  # These files are expected to exist in the package's inst/ directory
  # For standalone usage, adjust paths as needed
  
  # Load marker genes (40-gene signature)
  marker_genes_path <- system.file("inst", "marker_genes.csv", 
                                   package = "MASLD_subtyping")
  if (marker_genes_path == "") {
    # Fallback: look in current directory
    marker_genes_path <- "inst/marker_genes.csv"
  }
  
  marker_data <- read.csv(marker_genes_path, row.names = 1)
  marker_genes <- marker_data$gene_symbol
  
  message("Loaded ", length(marker_genes), " marker genes")
  
  # Load template for NTP method
  template_path <- system.file("inst", "ntp_templates.csv", 
                               package = "MASLD_subtyping")
  if (template_path == "") {
    template_path <- "inst/ntp_templates.csv"
  }
  templates <- read.csv(template_path, row.names = 1)
  
  # ---- 3. Subset Expression to Marker Genes ----
  # Identify which marker genes are present in the expression matrix
  genes_present <- intersect(marker_genes, rownames(expr_matrix))
  
  if (length(genes_present) == 0) {
    stop("None of the 40 marker genes found in the expression matrix. ",
         "Please ensure row names are HGNC gene symbols.")
  }
  
  if (length(genes_present) < 40) {
    warning("Only ", length(genes_present), " out of 40 marker genes present. ",
            "Classification accuracy may be reduced.")
  }
  
  # Subset expression matrix to marker genes
  expr_subset <- expr_matrix[genes_present, , drop = FALSE]
  
  # ---- 4. Perform Classification ----
  if (method == "NTP") {
    message("Performing classification using Nearest Template Prediction (NTP)...")
    
    results <- runNTP_standalone(
      expr       = as.data.frame(expr_subset),
      templates  = templates,
      scaleFlag  = TRUE,
      centerFlag = TRUE,
      doPlot     = do_plot,
      fig.name   = "MASLD_Subtype_Classification"
    )
    
  } else if (method == "PAM") {
    message("Performing classification using Prediction Analysis for Microarrays (PAM)...")
    
    # Load pre-trained PAM model
    pam_model_path <- system.file("inst", "pam_model.rds", 
                                  package = "MASLD_subtyping")
    if (pam_model_path == "") {
      pam_model_path <- "inst/pam_model.rds"
    }
    
    pam_model <- readRDS(pam_model_path)
    
    results <- runPAM_standalone(
      test.expr   = as.data.frame(expr_subset),
      train.model = pam_model,
      moic.res    = NULL
    )
  }
  
  # ---- 5. Prepare Output ----
  cluster_results <- results$clust.res
  
  # Ensure subtype labels are in CS1-CS4 format
  cluster_results$clust <- paste0("CS", cluster_results$clust)
  
  # Add sample IDs if not present
  if (!"sample_id" %in% colnames(cluster_results)) {
    cluster_results$sample_id <- rownames(cluster_results)
  }
  
  # ---- 6. Save Output Files ----
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    # Save cluster results
    write.csv(cluster_results,
              file.path(output_dir, "subtype_classification_results.csv"))
    
    # Save classification summary
    summary_df <- as.data.frame(table(cluster_results$clust))
    colnames(summary_df) <- c("Subtype", "Count")
    write.csv(summary_df,
              file.path(output_dir, "subtype_distribution_summary.csv"))
    
    message("Results saved to: ", output_dir)
  }
  
  # ---- 7. Return Results ----
  return(list(
    cluster_results = cluster_results,
    heatmap = if (do_plot) results$heatmap else NULL,
    summary = table(cluster_results$clust)
  ))
}


#' Standalone Implementation of NTP for MASLD Subtyping
#'
#' Internal function implementing Nearest Template Prediction for MASLD
#' subtype classification.
#'
#' @keywords internal
runNTP_standalone <- function(expr,
                              templates,
                              scaleFlag = TRUE,
                              centerFlag = TRUE,
                              doPlot = TRUE,
                              fig.name = "NTP_Heatmap") {
  
  # Check if required package is installed
  if (!requireNamespace("MOVICS", quietly = TRUE)) {
    warning("MOVICS package not installed. Falling back to correlation-based NTP.")
    return(runNTP_correlation(expr, templates, scaleFlag, centerFlag, doPlot))
  }
  
  # Use MOVICS package if available
  ntp_result <- MOVICS::runNTP(
    expr       = expr,
    templates  = templates,
    scaleFlag  = scaleFlag,
    centerFlag = centerFlag,
    doPlot     = doPlot,
    fig.name   = fig.name
  )
  
  return(ntp_result)
}


#' Alternative NTP Implementation Using Correlation
#'
#' @keywords internal
runNTP_correlation <- function(expr, templates, scaleFlag, centerFlag, doPlot) {
  
  # This is a simplified implementation for users without MOVICS
  # Core idea: compute correlation between each sample and each template,
  # assign subtype based on highest correlation
  
  # Scale and center if requested
  if (scaleFlag) {
    expr <- t(scale(t(expr)))
  }
  if (centerFlag) {
    expr <- expr - rowMeans(expr, na.rm = TRUE)
  }
  
  # Ensure genes match between expression and templates
  common_genes <- intersect(rownames(expr), rownames(templates))
  if (length(common_genes) == 0) {
    stop("No common genes between expression matrix and templates")
  }
  
  expr_matched <- expr[common_genes, ]
  templates_matched <- templates[common_genes, ]
  
  # Compute correlations
  cor_matrix <- cor(expr_matched, templates_matched, method = "pearson")
  
  # Assign each sample to the subtype with highest correlation
  best_subtype <- apply(cor_matrix, 1, function(x) {
    colnames(cor_matrix)[which.max(x)]
  })
  
  # Extract numeric subtype labels
  subtype_numeric <- as.numeric(gsub("CS", "", best_subtype))
  
  # Prepare results
  cluster_results <- data.frame(
    sample_id = rownames(cor_matrix),
    clust = subtype_numeric,
    confidence = apply(cor_matrix, 1, max),
    stringsAsFactors = FALSE
  )
  
  results <- list(
    clust.res = cluster_results,
    cor.matrix = cor_matrix
  )
  
  # Generate heatmap if requested
  if (doPlot) {
    if (requireNamespace("pheatmap", quietly = TRUE)) {
      results$heatmap <- pheatmap::pheatmap(
        t(cor_matrix),
        main = "NTP Classification Results",
        show_colnames = FALSE,
        cluster_cols = TRUE
      )
    } else {
      message("pheatmap package not installed; heatmap not generated")
    }
  }
  
  class(results) <- "NTPResult"
  return(results)
}


#' Standalone Implementation of PAM for MASLD Subtyping
#'
#' @keywords internal
runPAM_standalone <- function(test.expr, train.model, moic.res) {
  
  if (!requireNamespace("MOVICS", quietly = TRUE)) {
    stop("MOVICS package required for PAM classification. ",
         "Please install: install.packages('MOVICS')")
  }
  
  pam_result <- MOVICS::runPAM(
    train.expr = train.model$train.expr,
    moic.res   = moic.res,
    test.expr  = test.expr
  )
  
  return(pam_result)
}


# ---- Example Usage ----
#' Example: Classify samples from GSE193066
#'
#' @examples
#' \dontrun{
#' # Load expression data
#' exp <- read.csv("data/GSE193066_exp.csv", row.names = 1)
#' 
#' # Classify using NTP method
#' results <- classify_masld_subtype(
#'   expr_matrix = exp,
#'   method = "NTP",
#'   output_dir = "results/GSE193066",
#'   do_plot = TRUE
#' )
#' 
#' # View results
#' print(results$summary)
#' head(results$cluster_results)
#' }
NULL
