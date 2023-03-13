# StarTracer
## INTRODUCTION
A **beta-version**. This package will be able to used to search marker genes in an efficient way. It takes an `SeuratObject` or the direct output `data.frame` of `Seurat::FindAllMarkers()` function to search the ideal top N marker gene with a higher specificity level and efficiency.

The function also provides a intuitive way to descibe the specificity of the selected top N marker genes.
# INTSALLATION

starTracer could now be downloaded from github via the following command:

```
if (!requireNamespace("devtools", quietly = TRUE)){
install.packages("devtools")}

devtools::install_github("JerryZhang-1222/starTracer")
```

# USAGE

## Introduction

`starTracer` could be used to search Marker Genes from the single-cell/nuclear sequencing data. Marker genes are defined as the most up-regulated and most-specific genes in each cluster. `starTracer` will be able to generate the topN marker genes for each cluster in a effective way. The two most commonly used function are

1. `starTracer::searchMarker()`: the function takes in a `Seurat` object and calculate the marker genes for each cluster.
2. `starTracer::filterMarker()`: the function takes in a `data.frame` from the `Seurat::FindAllMArkers()` function, providing a new data.frame with an extra column named "order_trace", which could be used to reordering the marker genes from the `data.frame` the user provided.




