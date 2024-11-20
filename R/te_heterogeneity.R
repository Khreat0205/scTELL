#' Analyze TE heterogeneity using singleCellHaystack
#'
#' @param archr_proj ArchRProject object
#' @param matrix_name character Name of the TE matrix (default: "TEMatrix")
#' @param reduction character Name of dimensionality reduction (default: "UMAP")
#' @param grid_points numeric Number of grid points (default: 50)
#' @param p_threshold numeric P-value threshold for significance (default: 1e-10)
#' @return data.frame Significant heterogeneity results from Haystack analysis
#' @export
#'
#' @examples
#' \dontrun{
#' haystack_results <- analyze_te_heterogeneity(
#'   archr_proj = proj,
#'   reduction = "UMAP",
#'   grid_points = 50
#' )
#' }
analyze_te_heterogeneity <- function(
    archr_proj,
    matrix_name = "TEMatrix",
    reduction = "UMAP",
    grid_points = 50,
    p_threshold = 1e-10
) {
  # Input validation
  if (!methods::is(archr_proj, "ArchRProject")) {
    stop("archr_proj must be an ArchRProject object")
  }

  # Get TE matrix
  matrix <- ArchR::getMatrixFromProject(
    ArchRProj = archr_proj,
    useMatrix = matrix_name
  )
  matrix_data <- matrix@assays@data[[matrix_name]]
  rownames(matrix_data) <- matrix@elementMetadata$name

  # Get embedding
  embedding <- archr_proj@embeddings[[reduction]]$df

  # Match cell order
  matrix_data <- matrix_data[, match(rownames(embedding), colnames(matrix_data))]

  # Run haystack analysis
  res <- singleCellHaystack::haystack(
    embedding,
    matrix_data,
    grid.points = grid_points
  )

  # Get significant results
  res_sig <- singleCellHaystack::show_result_haystack(
    res.haystack = res,
    p.value.threshold = p_threshold
  )

  return(res_sig)
}

#' Plot TE accessibility on dimensional reduction
#'
#' @param archr_proj ArchRProject object
#' @param matrix_name character Name of the matrix (default: "TEMatrix")
#' @param feature_name character Name of the TE feature to plot
#' @param reduction character Name of dimensionality reduction (default: "UMAP")
#' @param use_imputation logical Whether to use imputation (default: TRUE)
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' p <- plot_te_heterogeneity(
#'   archr_proj = proj,
#'   feature_name = "L1HS_chr1_123456_123789_+",
#'   reduction = "UMAP"
#' )
#' }

plot_te_heterogeneity <- function(
    archr_proj,
    matrix_name = "TEMatrix",
    feature_name,
    reduction = "UMAP",
    use_imputation = TRUE
) {
  # Input validation
  if (!methods::is(archr_proj, "ArchRProject")) {
    stop("archr_proj must be an ArchRProject object")
  }

  # Add imputation weights if needed
  if (use_imputation && is.null(ArchR::getImputeWeights(archr_proj))) {
    archr_proj <- ArchR::addImputeWeights(archr_proj)
  }

  # Create plo87t
  p <- ArchR::plotEmbedding(
    ArchRProj = archr_proj,
    colorBy = matrix_name,
    name = feature_name,
    embedding = reduction,
    imputeWeights = if(use_imputation) ArchR::getImputeWeights(archr_proj) else NULL,
    randomize = FALSE
  )

  return(p)
}
