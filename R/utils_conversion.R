#' Convert ArchR project to Signac object
#'
#' @param archr_proj ArchRProject object
#' @param ref_version character Genome reference version (default: "hg38")
#' @param fragments_dir character Directory containing fragment files
#' @param peak_matrix matrix Peak matrix from ArchR project
#' @param fragments_from_cellranger logical Whether fragments are from cellranger (default: FALSE)
#' @param fragments_extension character File extension for fragment files
#' @param annotation Object containing gene annotations
#' @return Seurat object with Signac assay
#' @export
archr_to_signac <- function(
    archr_proj,
    ref_version = c("hg38", "hg19"),
    fragments_dir,
    peak_matrix,
    fragments_from_cellranger = FALSE,
    fragments_extension = NULL,
    annotation
) {
  # Input validation
  if (!methods::is(archr_proj, "ArchRProject")) {
    stop("archr_proj must be an ArchRProject object")
  }
  ref_version <- match.arg(ref_version)
  
  # Get samples
  samples <- unique(archr_proj@cellColData$Sample)
  
  # Process fragments path
  if (fragments_from_cellranger) {
    output_dir <- "/outs/"
    fragments_paths <- sapply(samples, function(sample) {
      paste0(fragments_dir, sample, output_dir, "fragments.tsv.gz")
    })
  } else {
    if (is.null(fragments_extension)) {
      stop("fragments_extension must be provided when fragments_from_cellranger = FALSE")
    }
    fragments_paths <- sapply(samples, function(sample) {
      paste0(fragments_dir, sample, fragments_extension)
    })
  }
  
  # Create Signac assay
  peak_assay <- Signac::CreateChromatinAssay(
    counts = peak_matrix,
    fragments = fragments_paths[1], # assuming single sample for now
    ranges = archr_proj@peakSet,
    genome = ref_version,
    annotation = annotation
  )
  
  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = peak_assay,
    assay = "peaks",
    meta.data = as.data.frame(archr_proj@cellColData)
  )
  
  return(seurat_obj)
}

#' Add dimensionality reductions from ArchR to Signac
#'
#' @param archr_proj ArchRProject object
#' @param seurat_obj Seurat object
#' @param reductions character vector Names of reductions to add
#' @param reduction_key character Key for reduction names
#' @return Updated Seurat object
#' @export
add_archr_reductions <- function(
    archr_proj,
    seurat_obj,
    reductions = c("UMAP", "IterativeLSI"),
    reduction_key = "archr_"
) {
  for (red in reductions) {
    if (red %in% names(archr_proj@embeddings)) {
      embed_df <- archr_proj@embeddings[[red]]$df
      rownames(embed_df) <- colnames(seurat_obj)
      
      seurat_obj[[paste0(reduction_key, tolower(red))]] <- 
        Seurat::CreateDimReducObject(
          embeddings = as.matrix(embed_df),
          assay = "peaks"
        )
    }
  }
  
  return(seurat_obj)
}