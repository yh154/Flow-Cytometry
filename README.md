# Flow-Cytometry
Scripts for Cytek data batch correction and clustering
``` r
Rscript flowcyto_bcc.R -h
```
Packages involved: 
``` r
suppressPackageStartupMessages({
require(flowCore)
require(CATALYST)
require(harmony)
require(data.table)
require(caret)
require(doParallel)
require(cowplot)
require(dplyr)
require(ggplot2)
})
```
