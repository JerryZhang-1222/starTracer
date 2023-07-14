# StarTracer

# 1 INTRODUCTION

A **beta-version**. This package will be able to used to search marker genes in an efficient way. It takes an `SeuratObject` or the direct output `data.frame` of `Seurat::FindAllMarkers()` function to search the ideal top N marker gene with a higher specificity level and efficiency.

The function also provides a intuitive way to descibe the specificity of the selected top N marker genes.

# 2 INTSALLATION

starTracer could now be downloaded from github via the following command:

```{r}
if (!requireNamespace("devtools", quietly = TRUE)){
install.packages("devtools")}

devtools::install_github("JerryZhang-1222/starTracer")
```

# 3 USAGE

## 3.1 Introduction

`starTracer` could be used to search Marker Genes from the single-cell/nuclear sequencing data. Marker genes are defined as the most up-regulated and most-specific genes in each cluster. `starTracer` will be able to generate the topN marker genes for each cluster in a effective way. The two most commonly used function are:

1.  `starTracer::searchMarker()`: a de-novo pipeline. The function takes in a `Seurat` object/Sparse matrix + Annotation/Average Expression Matrix and then calculate the marker genes for each cluster.
2.  `starTracer::filterMarker()`:an in-conjunction pipeline. The function takes in a `data.frame` from the `Seurat::FindAllMArkers()` function, providing a new data.frame with, which could be used to reordering the marker genes from the `data.frame` the user provided.
3.  
4.  

## 3.2 seachMrker

searchMarker is designed to take multiple kinds of input data. Users may provide input data as the following 3 formats:

1.  A Seurat object, with annotations.
2.  A Sparse Matrix of the single cell expression matrix and a annotation matrix.
3.  An average expression matrix.
4.  

For more details:

### 3.2.1 A Seurat Object

The object should meet certain requirements:

1.  Annotations should be provided in the Seurat object.
2.  High Variable Features should be found before running starTracer, if users prefer using "HVG" as genes to use.
3.  Integrated data should be avoided. If users have integrated data from multiple samples, please do not overwrite the original "RNA" slot.

```{r}
res <- searchMarker(
  x = x, #the input seurat object
  thresh.1 = 0.5,
  thresh.2 = NULL,
  method = "del_MI", 
  num = 2,
  gene.use = NULL,
  meta.data = NULL,
  ident.use = NULL
)
```

### 3.2.2 A Sparse Matrix

For users who may have processed the single cell experiment data in softwares other than Seurat, users are suggested to output their expression data as a Sparse Matrix(dgCMatrix) with cells in columns and features in rows and populated with the normalized data. An annotation matrix is also required with rows as cells corresponded with columns in the expression data.

```{r}
res <- searchMarker(
  x = x, #the input dgCMatrix
  thresh.1 = 0.5,
  thresh.2 = NULL,
  method = "del_MI",
  num = 2,
  gene.use = NULL,
  meta.data = meta.data, # the annotation matrix
  ident.use = NULL #the ident to use in the annotation matrix
)
```

### 3.2.3 An Average Expression Matrix

For users who have an average expression matrix with columns as clusters and rows as features, you may directly input the data.

```{r}
res <- searchMarker(
  x = x, #the average expression matrix
  thresh.1 = 0.5,
  thresh.2 = NULL,
  method = "del_MI",
  num = 2,
  gene.use = NULL,
  meta.data = NULL, # the annotation matrix
  ident.use = NULL #the ident to use in the annotation matrix
)
```
