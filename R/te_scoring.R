#' Compute TE accessibility scores for an ArchR project
#'
#' @param archr_proj ArchRProject object
#' @param te_loci GRanges object containing filtered TE annotations
#' @param name character Name for the matrix in ArchR project (default: "TEMatrix")
#' @param upstream_params list containing min and max upstream distances
#' @param downstream_params list containing min and max downstream distances
#' @param tile_size numeric Size of tiles for signal aggregation
#' @param decay_rate numeric Parameter for exponential decay weighting
#' @param gene_model character Model for signal weighting ("exp" or "uniform")
#' @param ceiling numeric Maximum value for accessibility scores
#' @param force logical Whether to force recomputation if matrix exists
#' @return Updated ArchRProject object with TE accessibility scores
#' @export
#'
#' @examples
#' \dontrun{
#' scored_proj <- compute_te_scores(
#'   archr_proj = proj,
#'   te_loci = filtered_catalog,
#'   upstream_params = list(min = 100, max = 1000),
#'   downstream_params = list(min = 100, max = 1000)
#' )
#' }
compute_te_scores <- function(
    archr_proj,
    te_loci,
    name = "TEMatrix",
    extend_upstream = 500,
    upstream_params = list(min = 100, max = 1000),
    downstream_params = list(min = 100, max = 1000),
    tile_size = 50,
    decay_rate = 1000,
    gene_model = c("exp", "uniform"),
    ceiling = 4,
    force = FALSE
) {

  # Input validation
  if (!is(archr_proj, "ArchRProject")) {
    stop("archr_proj must be an ArchRProject object")
  }
  if (!is(te_loci, "GRanges")) {
    stop("te_loci must be a GRanges object")
  }

  gene_model <- match.arg(gene_model)

  # Check if matrix already exists
  if (name %in% ArchR::getAvailableMatrices(archr_proj) && !force) {
    message("Matrix '", name, "' already exists. Use force = TRUE to recompute.")
    return(archr_proj)
  }

  # Create weighting function based on model
  weight_func <- if (gene_model == "exp") {
    sprintf('exp(-abs(x)/%s) + exp(-1)',decay_rate)
  } else {
    'abs(x)'
    # function(x) rep(1, length(x))
  }

  # Add TE score matrix
  archr_proj <- ArchR::addGeneScoreMatrix(
    input = archr_proj,
    genes = te_loci,
    matrixName = name,
    geneUpstream = extend_upstream,
    extendUpstream = c(upstream_params$min, upstream_params$max),
    extendDownstream = c(downstream_params$min, downstream_params$max),
    geneModel = weight_func,
    tileSize = tile_size,
    ceiling = ceiling,
    force = force
  )

  return(archr_proj)
}

#' Extract TE accessibility matrix from ArchR project
#'
#' @param archr_proj ArchRProject object
#' @param matrix_name character Name of the matrix to extract
#' @param scale logical Whether to scale the matrix
#' @param scale_factor numeric Factor for scaling
#' @return SummarizedExperiment object containing TE accessibility scores
#' @export
#'
#' @examples
#' \dontrun{
#' te_mat <- extract_te_matrix(
#'   archr_proj = scored_proj,
#'   matrix_name = "TEMatrix"
#' )
#' }
extract_te_matrix <- function(
    archr_proj,
    matrix_name = "TEMatrix",
    scale = FALSE,
    scale_factor = 10000
) {

  if (!matrix_name %in% ArchR::getAvailableMatrices(archr_proj)) {
    stop("Matrix '", matrix_name, "' not found in ArchR project")
  }

  # Extract matrix
  te_mat <- ArchR::getMatrixFromProject(
    ArchRProj = archr_proj,
    useMatrix = matrix_name
  )

  # Scale if requested
  if (scale) {
    assay_data <- SummarizedExperiment::assay(te_mat)
    lib_sizes <- colSums(assay_data)
    assay_data <- t(t(assay_data) / lib_sizes * scale_factor)
    assay(te_mat) <- assay_data
  }

  return(te_mat)
}

#' Filter cells based on TE accessibility
#'
#' @param archr_proj ArchRProject object
#' @param matrix_name character Name of the TE matrix
#' @param min_counts numeric Minimum number of counts required
#' @param min_detected numeric Minimum number of TEs detected
#' @return Vector of cell names passing filters
#' @export
#'
#' @examples
#' \dontrun{
#' valid_cells <- filter_te_cells(
#'   archr_proj = scored_proj,
#'   min_counts = 100,
#'   min_detected = 10
#' )
#' }
filter_te_cells <- function(
    archr_proj,
    matrix_name = "TEMatrix",
    min_counts = NULL,
    min_detected = NULL
) {

  te_mat <- ArchR::getMatrixFromProject(
    ArchRProj = archr_proj,
    useMatrix = matrix_name
  )

  assay_data <- SummarizedExperiment::assay(te_mat)

  # Calculate metrics
  total_counts <- colSums(assay_data)
  detected_tes <- colSums(assay_data > 0)

  # Apply filters
  keep_cells <- rep(TRUE, ncol(assay_data))

  if (!is.null(min_counts)) {
    keep_cells <- keep_cells & (total_counts >= min_counts)
  }

  if (!is.null(min_detected)) {
    keep_cells <- keep_cells & (detected_tes >= min_detected)
  }

  valid_cells <- colnames(assay_data)[keep_cells]

  # Return statistics
  stats <- data.frame(
    total_cells = ncol(assay_data),
    passing_cells = sum(keep_cells),
    mean_counts = mean(total_counts),
    mean_detected = mean(detected_tes)
  )

  message("Filtering statistics:\n")
  print(stats)

  return(valid_cells)
}

#' Calculate element-level accessibility scores
#'
#' @param te_matrix SummarizedExperiment TE accessibility matrix
#' @param te_granges GRanges object containing TE annotations
#' @return data.frame with element-level accessibility scores for each cell
#' @export
#'
#' @examples
#' \dontrun{
#' element_scores <- calculate_element_scores(
#'   te_matrix = te_mat,
#'   te_granges = te_loci
#' )
#' }
calculate_element_scores <- function(
    te_matrix,
    te_granges
) {
  # Input validation
  if (!methods::is(te_matrix, "SummarizedExperiment")) {
    stop("te_matrix must be a SummarizedExperiment object")
  }
  if (!methods::is(te_granges, "GRanges")) {
    stop("te_granges must be a GRanges object")
  }

  # Get accessibility matrix and transpose
  scaled_mat <- t(SummarizedExperiment::assay(te_matrix))
  colnames(scaled_mat) <- SummarizedExperiment::rowData(te_matrix)$name

  # Get unique element names
  element_names <- unique(gsub("_chr.*", "", te_granges$symbol))

  # Initialize results list
  element_scores <- list()

  # Calculate scores for each element type
  for (element in element_names) {
    # Get loci for current element
    locus_names <- te_granges$symbol[grep(sprintf("%s_chr", element),
                                          te_granges$gene_id)]

    if (length(locus_names) == 0) next

    if (length(locus_names) == 1) {
      element_scores[[element]] <- scaled_mat[, colnames(scaled_mat) %in% locus_names,
                                              drop = FALSE]
    } else {
      element_scores[[element]] <- rowSums(scaled_mat[, colnames(scaled_mat) %in% locus_names,
                                                      drop = FALSE],
                                           na.rm = TRUE)
    }
  }

  # Convert to data.frame
  score_df <- as.data.frame(do.call(cbind, element_scores))
  colnames(score_df) <- paste0("score_", gsub("-", "_", element_names))
  rownames(score_df) <- rownames(scaled_mat)

  return(score_df)
}

#' Add element scores to ArchR project
#'
#' @param archr_proj ArchRProject object
#' @param score_df data.frame containing element scores
#' @param overwrite logical whether to overwrite existing scores
#' @return Updated ArchRProject object
#' @export
#'
#' @examples
#' \dontrun{
#' updated_proj <- add_element_scores(
#'   archr_proj = proj,
#'   score_df = element_scores
#' )
#' }
add_element_scores <- function(
    archr_proj,
    score_df,
    overwrite = FALSE
) {

  if (!methods::is(archr_proj, "ArchRProject")) {
    stop("archr_proj must be an ArchRProject object")
  }

  # Check cell names match
  if (!all(rownames(score_df) == ArchR::getCellNames(archr_proj))) {
    stop("Cell names in score_df do not match ArchR project")
  }

  # Add scores to cellColData
  for (score_name in colnames(score_df)) {
    if (score_name %in% names(archr_proj@cellColData) && !overwrite) {
      message(score_name, " already exists. Set overwrite = TRUE to replace.")
      next
    }
    archr_proj <- ArchR::addCellColData(
      ArchRProj = archr_proj,
      data = score_df[[score_name]],
      name = score_name,
      cells = rownames(score_df)
    )
  }

  return(archr_proj)
}
