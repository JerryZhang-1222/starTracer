---
title: "starTracer-vignette"
output: rmarkdown::html_vignette
date: "`r Sys.Date()`"
author: 
- name: Feiyang Zhang
  email: fyzhang.neu@gmail.com
  affiliation: Wuhan University
vignette: >
  %\VignetteIndexEntry{starTracer-vignette}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Overview

starTracer is an R package for marker gene identification in single-cell RNA-seq data analysis. This package comprises two primary functional modules: "searchMarker" and "filterMarker". The "searchMarker" function operates independently, taking an expression matrix as input and generating a marker gene matrix as output. On the other hand, the "filterMarker" function is tailored to work as a complementary pipeline to the Seurat "FindAllMarkers" function, offering a more accurate list of marker genes for each cluster in conjunction with Seurat results.

# Testing Scripts of searchMarker

Please refer to the following links for details in our original research article. You may also find the basic usage examples here:

## Speed, Specificity and Accuracy

In the following 3 scripts, we test the specificity, calculation speed and accuracy in three samples:

[human brain](https://jerryzhang-1222.github.io/brain.html)

[human heart](https://jerryzhang-1222.github.io/heart_big.html)

[mouse kidney](https://jerryzhang-1222.github.io/kidney.html)

## Ability to Find Marker Genes in Different Cluster Levels

We here use the human brain sample, with the annotation of "bio_clust", "major_clust" and "sub_clust", to test `starTracer::searchMarker`'s ability to find maker genes across different annotation levels.

[different cluster level](https://jerryzhang-1222.github.io/different_cluster_level.html)

## Parameter Control

Here we show the different results of using different thresh.2 (denoted as S2 in the original research article). You may refer to this results to better understand the influence of thresh.2 on your results.

[different threshold](https://jerryzhang-1222.github.io/Different_thresh.html)

# Testing Scripts of filterMarker

Here we show the results of filterMarker. filterMarker takes the output matrix/data.frame of Seurat's FindAllMarkers function. Note that calculating with samples of big cell numbers will be time consumable.

[filterMarker test](https://jerryzhang-1222.github.io/FilterMarkerTest_v3.html)
