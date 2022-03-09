
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scUtils

<!-- badges: start -->
<!-- badges: end -->

Single-cell utility functions by Lukas Steuernagel

## Installation

Install scUtils using:

``` r
#devtools::install_github("lsteuernagel/scUtils")
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
#>  [1] "apply_DoubletFinder"           "determine_cluster_resolution" 
#>  [3] "downsample_balanced_iterative" "gene_pct_cluster"             
#>  [5] "identify_variable_features"    "infer_sex"                    
#>  [7] "match_sample_names"            "rasterize_ggplot"             
#>  [9] "Read10xFormat"                 "ReadDGEFormat"                
#> [11] "seurat_recipe"                 "writeList_to_JSON"
```

Check the man pages for details!

## Example

Still need to add examples.
