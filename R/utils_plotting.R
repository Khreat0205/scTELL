#' Plot TE coverage using Signac
#'
#' @param seurat_obj Seurat object
#' @param te_regions GRanges object containing TE regions
#' @param group_by character Column name for grouping
#' @param extend_upstream numeric Bases to extend upstream
#' @param extend_downstream numeric Bases to extend downstream
#' @return ggplot object
#' @export
plot_te_coverage <- function(
    seurat_obj,
    te_regions,
    group_by = NULL,
    extend_upstream = 3000,
    extend_downstream = 3000
) {
  # Create coverage plot
  p_coverage <- Signac::CoveragePlot(
    object = seurat_obj,
    region = te_regions,
    group.by = group_by,
    extend.upstream = extend_upstream,
    extend.downstream = extend_downstream
  )
  
  # Add annotation plot
  p_annotation <- Signac::AnnotationPlot(
    object = seurat_obj,
    region = te_regions
  )
  
  # Combine plots
  p <- patchwork::wrap_plots(
    p_coverage, p_annotation,
    ncol = 1,
    heights = c(3, 1)
  )
  
  return(p)
}

#' Create heatmap of TE accessibility
#'
#' @param te_matrix matrix TE accessibility matrix
#' @param group_info factor Group labels for columns
#' @param scale logical Whether to scale data
#' @param show_rownames logical Whether to show row names
#' @return ComplexHeatmap object
#' @export
plot_te_heatmap <- function(
    te_matrix,
    group_info = NULL,
    scale = TRUE,
    show_rownames = TRUE
) {
  # Scale data if requested
  if (scale) {
    mat_scaled <- t(scale(t(te_matrix)))
  } else {
    mat_scaled <- te_matrix
  }
  
  # Create annotation if group info provided
  if (!is.null(group_info)) {
    ha <- ComplexHeatmap::HeatmapAnnotation(
      group = group_info
    )
  } else {
    ha <- NULL
  }
  
  # Create heatmap
  ht <- ComplexHeatmap::Heatmap(
    matrix = mat_scaled,
    name = "Accessibility",
    show_row_names = show_rownames,
    top_annotation = ha,
    cluster_columns = TRUE,
    cluster_rows = TRUE
  )
  
  return(ht)
}

#' Plot QC metrics for TE accessibility
#'
#' @param archr_proj ArchRProject object
#' @param matrix_name character Name of TE matrix
#' @return List of ggplot objects
#' @export
plot_te_qc <- function(
    archr_proj,
    matrix_name = "TEMatrix"
) {
  # Get matrix
  te_mat <- ArchR::getMatrixFromProject(
    ArchRProj = archr_proj,
    useMatrix = matrix_name
  )
  
  # Calculate QC metrics
  total_counts <- Matrix::colSums(assay(te_mat))
  detected_tes <- Matrix::colSums(assay(te_mat) > 0)
  
  # Create plots
  p1 <- ggplot2::ggplot(
    data.frame(counts = total_counts),
    ggplot2::aes(x = counts)
  ) +
    ggplot2::geom_histogram() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Distribution of TE counts per cell")
  
  p2 <- ggplot2::ggplot(
    data.frame(detected = detected_tes),
    ggplot2::aes(x = detected)
  ) +
    ggplot2::geom_histogram() +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "Distribution of detected TEs per cell")
  
  return(list(counts = p1, detected = p2))
}