# StarTracer

[![DOI](https://zenodo.org/badge/611591370.svg)](https://zenodo.org/doi/10.5281/zenodo.13364966)

# 1 INTRODUCTION

This package will be able to used to search marker genes in an efficient way, which will accelerate your marker gene search process hundreds to thousands times faster.The function also provides a intuitive way to descibe the specificity of the selected top N marker genes.


starTracer takes an `SeuratObject` as input or the direct output `data.frame` of `Seurat::FindAllMarkers()` function to search the ideal top N marker gene with a higher specificity level and efficiency.

Please cite `doi: https://doi.org/10.1101/2023.09.21.558919` if starTracer has offered your help during your marker gene search process.

## Instruction of How to Choose the Functions in starTracer:
1. If you have a Seurat Object: you may directly put it into function `searchMarker`, please see details in usage or vignette.
2. If you have any trouble with running with your Seurat Object, you may a). Utilizing `AverageExpress()`  from Seurat to generate a pseudo-bulk expression matrix, then impot the matrix into searchMarker function. b) Output your dgCMatrix, then import your matrix into `searchMarker` function.
3. If you already have your Seurat's output, you can directly import that into `filterMarker` function to get a refined result.
   Note that assuming using a same expression data, results from 1 and 2 will be the same while result 3 will differ.

## Special Instruction for Python or other none-R users:
We suggest you output your data as an average expression matrix (Pseudo-bulk Expression Matrix), then import it into R.
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

## 3.2 searchMrker

searchMarker is designed to take multiple kinds of input data. Users may provide input data as the following 3 formats:

1.  A Seurat object, with annotations.
2.  A Sparse Matrix of the single cell expression matrix and a annotation matrix.
3.  An average expression matrix.

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
  method = "pos", 
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
  method = "pos",
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
  method = "pos",
  num = 2,
  gene.use = NULL,
  meta.data = NULL, # the annotation matrix
  ident.use = NULL #the ident to use in the annotation matrix
)
```
## 3.3 filterMarker, Automatically Refine Your Marker Gene Results from Seurat.
Preparing your FIndAllMarkers output as a data.frame, the Seurat Object and run the following sripts, which might take a while before getting the refined result.
```{r}
res.filt <- starTracer::filterMarker(x = <FindAllMarkers output>,
                                     ident.use = <your idents>,
                                     <your Seurat Object>,
                                     num = 5)
```
An additional column will be added as `pct.pos`, after which, you may
1. Filter out marker genes according to a specific p-value and Fold-change.
2. For each cluster, sort the marker gene by descending `pct.pos`.
   
## 3.4 NEW FUNCTION (Alpha Version): Find Negative Markers

We now offer users a solution to find negative markers. Imply by specifying method = "neg"
```{r}
res <- searchMarker(
  x = x,
  thresh.1 = 0.5,
  thresh.2 = NULL,
  method = "neg",
  num = 2,
  gene.use = NULL,
  meta.data = NULL, # the annotation matrix
  ident.use = NULL #the ident to use in the annotation matrix
)
```
