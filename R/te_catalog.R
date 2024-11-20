#' Create TE catalog from RepeatMasker table
#'
#' @param rmsk_table data.frame containing RepeatMasker annotations
#' @param element_types character vector of TE types to include (e.g., "L1HS", "SVA_A")
#' @param min_length numeric minimum length requirement for TEs
#' @param genome character genome build version ("hg38" or "hg19")
#' @param remove_alt_chroms logical whether to remove alternative chromosomes (default TRUE)
#' @return GRanges object containing filtered TE annotations
#' @export
#'
#' @examples
#' \dontrun{
#' rmsk <- read.table("rmsk.txt")
#' l1_catalog <- create_te_catalog(
#'   rmsk_table = rmsk,
#'   element_types = c("L1HS", "L1PA2", "L1PA3"),
#'   min_length = 5000
#' )
#' }
create_te_catalog <- function(rmsk_table,
                              element_types = NULL,
                              min_length = NULL,
                              genome = c("hg38", "hg19"),
                              remove_alt_chroms = TRUE) {

  # Input validation
  if(!is.data.frame(rmsk_table)) {
    stop("rmsk_table must be a data.frame")
  }

  genome <- match.arg(genome)

  # Filter by element types if specified
  if(!is.null(element_types)) {
    rmsk_table <- rmsk_table[rmsk_table$repName %in% element_types,]
    if(nrow(rmsk_table) == 0) {
      stop("No elements found matching specified element_types")
    }
  }

  # Filter by length if specified
  if(!is.null(min_length)) {
    element_lengths <- abs(rmsk_table$genoEnd - rmsk_table$genoStart)
    rmsk_table <- rmsk_table[element_lengths >= min_length,]
    if(nrow(rmsk_table) == 0) {
      stop("No elements found meeting minimum length requirement")
    }
  }

  # Remove alternative chromosomes if specified
  if(remove_alt_chroms) {
    rmsk_table <- rmsk_table[!grepl("_", rmsk_table$genoName),]
  }

  # Create GRanges object
  te_granges <- GenomicRanges::GRanges(
    seqnames = rmsk_table$genoName,
    ranges = IRanges::IRanges(
      start = rmsk_table$genoStart,
      end = rmsk_table$genoEnd
    ),
    strand = rmsk_table$strand,
    gene_id = paste(rmsk_table$repName,
                    rmsk_table$genoName,
                    rmsk_table$genoStart,
                    rmsk_table$genoEnd,
                    rmsk_table$strand,
                    sep = "_"),
    symbol = paste(rmsk_table$repName,
                   rmsk_table$genoName,
                   rmsk_table$genoStart,
                   rmsk_table$genoEnd,
                   rmsk_table$strand,
                   sep = "_")
  )

  return(te_granges)
}

#' Combine multiple TE catalogs
#'
#' @param catalog_list list of GRanges objects containing TE annotations
#' @return Combined GRanges object
#' @export
#'
#' @examples
#' \dontrun{
#' combined_catalog <- combine_te_catalogs(
#'   list(l1_catalog, sva_catalog)
#' )
#' }
combine_te_catalogs <- function(catalog_list) {
  # Input validation
  if(!is.list(catalog_list)) {
    stop("catalog_list must be a list")
  }

  if(!all(sapply(catalog_list, is, "GRanges"))) {
    stop("All elements in catalog_list must be GRanges objects")
  }

  # Combine catalogs
  combined <- do.call(c, catalog_list)

  # Remove any potential duplicates
  combined <- unique(combined)

  return(combined)
}

#' Filter TE loci based on peaks and blacklist
#'
#' @param te_catalog GRanges object containing TE annotations
#' @param peak_regions GRanges object containing ATAC-seq peaks
#' @param blacklist GRanges object containing blacklist regions
#' @param upstream_extension numeric bases to extend upstream for regulatory regions
#' @return GRanges object containing filtered TE annotations
#' @export
#'
#' @examples
#' \dontrun{
#' filtered_catalog <- filter_te_loci(
#'   te_catalog = combined_catalog,
#'   peak_regions = peaks,
#'   blacklist = blacklist_regions,
#'   upstream_extension = 1000
#' )
#' }
filter_te_loci <- function(te_catalog,
                           peak_regions = NULL,
                           blacklist = NULL,
                           upstream_extension = 1000) {

  # Input validation
  if(!is(te_catalog, "GRanges")) {
    stop("te_catalog must be a GRanges object")
  }

  if(!is.null(peak_regions) && !is(peak_regions, "GRanges")) {
    stop("peak_regions must be a GRanges object")
  }

  if(!is.null(blacklist) && !is(blacklist, "GRanges")) {
    stop("blacklist must be a GRanges object")
  }

  # Create extended regions for regulatory analysis
  te_extended <- GenomicRanges::GRanges(
    seqnames = seqnames(te_catalog),
    ranges = IRanges::IRanges(
      start = ifelse(strand(te_catalog) == "-",
                     start(te_catalog),
                     start(te_catalog) - upstream_extension),
      end = ifelse(strand(te_catalog) == "-",
                   end(te_catalog) + upstream_extension,
                   end(te_catalog))
    ),
    strand = strand(te_catalog)
  )

  # Copy metadata
  mcols(te_extended) <- mcols(te_catalog)

  # Filter based on peak overlap
  if(!is.null(peak_regions)) {
    peak_overlap <- GenomicRanges::findOverlaps(te_extended, peak_regions)
    te_extended <- te_extended[unique(queryHits(peak_overlap))]
    te_catalog <- te_catalog[unique(queryHits(peak_overlap))]
  }

  # Remove blacklisted regions
  if(!is.null(blacklist)) {
    blacklist_overlap <- GenomicRanges::findOverlaps(te_extended, blacklist)
    if(length(blacklist_overlap) > 0) {
      te_catalog <- te_catalog[-unique(queryHits(blacklist_overlap))]
    }
  }

  return(te_catalog)
}
