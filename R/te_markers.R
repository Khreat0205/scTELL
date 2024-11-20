#' Find TE markers for specified groups
#'
#' @param archr_proj ArchRProject object
#' @param matrix_name character Name of the TE matrix
#' @param group_by character Column name in cellColData for grouping
#' @param test_method character Statistical test method ("wilcoxon", "t", "LR")
#' @param fdr_threshold numeric FDR threshold for significance
#' @param log2fc_threshold numeric Log2FC threshold for significance
#' @param min_pct numeric Minimum percentage of cells expressing TE
#' @param min_counts numeric Minimum counts for filtering cells (passed to filter_te_cells)
#' @param min_detected numeric Minimum detected TEs for filtering cells (passed to filter_te_cells)
#' @return List of marker statistics for each group
#' @export
find_te_markers <- function(
    archr_proj,
    matrix_name = "TEMatrix",
    group_by,
    test_method = c("wilcoxon", "ttest", "binomial"),
    fdr_threshold = 0.001,
    log2fc_threshold = 1,
    min_pct = 0.1,
    min_counts = NULL,
    min_detected = NULL
) {
  # Input validation
  if (!methods::is(archr_proj, "ArchRProject")) {
    stop("archr_proj must be an ArchRProject object")
  }

  test_method <- match.arg(test_method)

  # Filter cells using existing function
  valid_cells <- filter_te_cells(
    archr_proj = archr_proj,
    matrix_name = matrix_name,
    min_counts = min_counts,
    min_detected = min_detected
  )

  # Get markers using filtered cells
  markers <- ArchR::getMarkerFeatures(
    ArchRProj = ArchR::subsetCells(archr_proj, cellNames=valid_cells),
    useMatrix = matrix_name,
    groupBy = group_by,
    testMethod = test_method,
    bias = c("TSSEnrichment", "log10(nFrags)")
  )
  print(markers)
  # Get significant markers
  markers_sign <- ArchR::getMarkers(
    markers,
    cutOff = sprintf("FDR <= %s & Log2FC >= %s", fdr_threshold, log2fc_threshold)
  )

  return(markers_sign)
}

#' Validate markers across datasets
#'
#' @description
#' Validates markers by finding common elements between discovery and validation datasets,
#' using a cell type matching table to handle different naming conventions.
#'
#' @param discovery_markers List of markers from discovery dataset
#' @param validation_markers List of markers from validation dataset
#' @param matching_table data.frame with following columns:
#'   \itemize{
#'     \item cell_type: Cell type names in discovery dataset
#'     \item cell_type_val: Corresponding cell type names in validation dataset
#'   }
#' @return List of validated markers, where names correspond to discovery dataset cell types
#'
#' @examples
#' \dontrun{
#' matching_df <- data.frame(
#'   cell_type = c("B cells", "CD8 T", "NK"),
#'   cell_type_val = c("B", "CD8_T", "NK_cell")
#' )
#' validated <- validate_te_markers(
#'   discovery_markers = markers1,
#'   validation_markers = markers2,
#'   matching_table = matching_df
#' )
#' }
#' @export
validate_te_markers <- function(
    discovery_markers,
    validation_markers,
    matching_table
) {
  # Input validation
  if (!all(c("cell_type", "cell_type_val") %in% colnames(matching_table))) {
    stop("matching_table must contain 'cell_type' and 'cell_type_val' columns")
  }

  # Match cell type names
  matched_names <- matching_table$celltype_val[match(
    names(discovery_markers),
    matching_table$celltype
  )]

  # Initialize results list
  valid_markers <- list()

  # Find common markers
  for (i in seq_along(matched_names)) {
    if (is.na(matched_names[i])) next

    discovery_name <- names(discovery_markers)[i]
    validation_name <- matched_names[i]

    if (validation_name %in% names(validation_markers)) {
      common_markers <- intersect(
        discovery_markers[[discovery_name]]$name,
        validation_markers[[validation_name]]$name
      )

      if (length(common_markers) > 0) {
        valid_markers[[discovery_name]] <- common_markers
      }
    }
  }

  return(valid_markers)
}


