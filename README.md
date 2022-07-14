
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scUtils

<!-- badges: start -->
<!-- badges: end -->

Single-cell utility functions by Lukas Steuernagel

## Installation

Install scUtils using:

``` r
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")
remotes::install_github("lsteuernagel/scUtils")
```

Depends on Seurat, some tidyverse and other data management packages.

## Function overview

Load package:

``` r
library(scUtils)
```

List of functions:

``` r
sort(as.character(lsf.str("package:scUtils")))
#>  [1] "apply_DoubletFinder"           "CalculateMultScore"           
#>  [3] "conserved_correlations"        "determine_cluster_resolution" 
#>  [5] "downsample_balanced_iterative" "feature_statistics"           
#>  [7] "find_children"                 "gene_pct_cluster"             
#>  [9] "identify_variable_features"    "infer_sex"                    
#> [11] "match_sample_names"            "rasterize_ggplot"             
#> [13] "Read10xFormat"                 "ReadDGEFormat"                
#> [15] "seurat_recipe"                 "writeList_to_JSON"
```

Check the man pages for details!

## Examples

Still need to add examples…
