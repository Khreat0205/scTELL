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

# Define target TEs
te_regions <- defineTargetTEs(
  genome = "hg38",
  te_families = c("L1HS", "SVA", "HERVK"),
  min_length = c(5000, 700, 5000)
)

# Calculate TE accessibility scores
te_matrix <- calculateTEScores(
  project = proj,
  te_regions = te_regions,
  extend_upstream = 1000
)

# Perform cell type-specific analysis
te_markers <- findTEMarkers(
  te_matrix = te_matrix,
  group_by = "cell_type",
  test_method = "wilcoxon"
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
