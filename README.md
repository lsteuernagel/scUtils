
<!-- README.md is generated from README.Rmd. Please edit that file -->

# scUtils

<!-- badges: start -->
<!-- badges: end -->

Single-cell utility functions

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
#>  [7] "find_children"                 "FindDEG_nebula"               
#>  [9] "gene_pct_cluster"              "identify_variable_features"   
#> [11] "infer_sex"                     "match_sample_names"           
#> [13] "rasterize_ggplot"              "Read10xFormat"                
#> [15] "ReadDGEFormat"                 "run_nebula"                   
#> [17] "seurat_recipe"                 "writeList_to_JSON"
```

Check the man pages for details!
