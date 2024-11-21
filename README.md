# scTELL

## Overview
scTELL (single-cell Transposable Element Locus-Level analysis) is a computational framework designed to analyze transposable element (TE) accessibility patterns at locus-level resolution in single-cell ATAC-seq data. scTELL enables comprehensive exploration of TE activity across different cell populations and can identify cell type-specific and heterogeneous TE accessibility patterns.

## Features
- Locus-level quantification of TE accessibility in scATAC-seq data
- Cell type-specific TE analysis with built-in validation framework
- Cluster-free analysis for capturing tumor heterogeneity 
- Integration with bulk ATAC-seq data for validation
- Comprehensive visualization tools

## Installation

```R
# Install devtools if not already installed
if (!require("devtools")) {
  install.packages("devtools")
}

# Install scTELL
devtools::install_github("Khreat0205/scTELL")
```

## Dependencies
- R >= 4.1.0
- ArchR >= 1.0.1
- Signac >= 1.7.0
- singleCellHaystack >= 0.3.3
- GenomicRanges >= 1.46.1
- data.table >= 1.14.2

## Quick Start

```R
library(scTELL)

# Load your data
proj <- loadArchRProject("path/to/ArchR/project")

rmsk <- read.table("path/to/rmsk/table")

# Define target TEs
l1_catalog <- create_te_catalog(
  rmsk_table = rmsk,
  element_types = c("L1HS", "L1PA2", "L1PA3"),
  min_length = 5000
)

filtered_catalog <- filter_te_loci(
  te_catalog = combined_catalog,
  peak_regions = peaks,
  blacklist = blacklist_regions,
  upstream_extension = 1000
)

# Calculate TE accessibility scores
scored_proj <- compute_te_scores(
  archr_proj = proj,
  te_loci = filtered_catalog,
  upstream_params = list(min = 100, max = 1000),
  downstream_params = list(min = 100, max = 1000)
)

te_mat <- extract_te_matrix(
  archr_proj = scored_proj,
  matrix_name = "TEMatrix"
)

# Perform cell type-specific analysis
te_markers <- find_te_markers(
  archr_proj = proj,
  group_by = "cell_type",
  test_method = "wilcoxon"
)

# Perform cluster-free analysis
haystack_results = analyze_te_heterogeneity(
  archr_proj = proj,
  reduction = "UMAP",
  grid_points = 50
)


```


## Citation
If you use scTELL in your research, please cite:
```
[---]
```

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact
For questions and feedback:
- Submit an issue on GitHub
- Email: scientist0205@snu.ac.kr

## Contributing
We welcome contributions! Please see our [contributing guidelines](CONTRIBUTING.md) for details.
