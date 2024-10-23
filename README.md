# FSBLUP

<!-- badges: start -->

<!-- badges: end -->

## **F**usion **S**imilarity **M**atrix **B**est **L**inear **U**nbiased **P**rediction

## CONTENES

-   [OVERVIEW](#overview)
-   [GETTING STARTED](#getting-started)
    -   [Installation](#installation)
    -   [Test datasets](#test-datasets)
-   [Data input](#data-input)
    -   [Phenotype data](#phenotype-data)
    -   [Kinship data](#kinship-data)
    -   [Pedigree data](#pedigree-data)
    -   [Genotype data](#genotype-data)
    -   [Intermediate omics data](#intermediate-omics-data)
-   [USAGE](#usage)
    -   [Basic](#basic)
    -   [Advanced](#advanced)
-   [OUTPUT](#output)

------------------------------------------------------------------------

## OVERVIEW

`FSBLUP` is a novel strategy of fusion similarity matrix based on genomic and intermediate omics information, to estimate the unique genetic correlation of different data with target phenotype systematically and adapt the genetic architecture of complex traits jointly, by a machine learning-based method incorporating cross-validation, grid search, and adaptive bisection algorithms.

## GETTING STARTED

***`FSBLUP`*** is implemented in R, and for faster computation with large datasets, it is recommended to link R with either the Math Kernel Library (MKL) or OpenBLAS. Both MKL and OpenBLAS can accelerate the BLAS/LAPACK library using multi-threading capabilities, significantly reducing computational time. Integrating MKL or OpenBLAS with R can improve performance by automatically optimizing matrix operations across multiple threads. Detailed instructions are available on how to link R with these libraries to achieve enhanced computational efficiency.

> 1. [how to link MKL with R](https://www.intel.com/content/www/us/en/developer/articles/technical/quick-linking-intel-mkl-blas-lapack-to-r.html)

> 2. [how to link OpenBLAS with R](https://stackoverflow.com/questions/74086774/build-r-with-openblas)

### Installation

You can install the development version of `FSBLUP` like so:

``` r
#if "devtools" isn't installed, please "install.packages(devtools)" first.
devtools::install_github("CAU-TeamLiuJF/FSBLUP")
```

After installed successfully, the `FSBLUP` package can be loaded by typing

``` r
library(FSBLUP)
```

Typing `?FSBLUP` could get the details of all parameters.

### Test datasets

The example data can be downloaded by typing:

``` bash
wget https://github.com/CAU-TeamLiuJF/FSBLUP/releases/download/example/example.zip
unzip example.zip
```

## DATA INPUT

### Phenotype data

The phenotype should be a `data.frame` contain a header row, e.g. trait names. The missing values should be denoted by NA, which will be treated as candidates. Notice that only the numeric values are allowed and the characters will not be recognized. The first column should be IDs.

> phenotype

|    ID     |     A     |  B  | T1  |    T2    |
|:---------:|:---------:|:---:|:---:|:--------:|
| 0.224992  | 0.224991  |  1  | NA  | 0.285427 |
| -0.974543 | -0.974542 | NA  |  0  | 0.195909 |
|    NA     |    NA     | NA  | NA  |    NA    |
|    ...    |    ...    | ... | ... |   ...    |

### Kinship data

`FSBLUP` accept kinship matrix instead of original data as input. The matrices should be ${n}\times{n}$ that represents the relationship among individuals.

It could be also supplied by the users, however, in this case, the order of individuals in either row or column should be the same as phenotype file, the column and row names should be supplied.

> kinship

|     |  1   |   2   |   3   |   4   |   5   |   6   | ... |
|:---:|:----:|:-----:|:-----:|:-----:|:-----:|:-----:|:---:|
|  1  | 0.90 | 0.17  | 0.17  | 0.02  | 0.17  | 0.19  | ... |
|  2  | 0.17 | 1.59  | 0.99  | 0.17  | -0.10 | -0.21 | ... |
|  3  | 0.17 | 0.99  | 1.05  | 0.11  | -0.08 | 0.00  | ... |
|  4  | 0.02 | 0.17  | 0.11  | 0.60  | -0.01 | -0.13 | ... |
|  5  | 0.17 | -0.10 | -0.08 | -0.01 | 0.94  | 0.08  | ... |
|  6  | 0.19 | -0.21 | 0.00  | -0.13 | 0.08  | 0.96  | ... |
| ... | ...  |  ...  |  ...  |  ...  |  ...  |  ...  | ... |

### Pedigree data

The pedigree should be a `data.frame` or `matrix` structure, which will be used to determine additive genetic relationships. The pedigree should include three columns:

1.  ID
    -   The individual number for random effects (genetic effects) in the model.
2.  Sire ID
    -   The sire number; if unknown, use 0.
3.  Dam ID
    -   The dam number; if unknown, use 0.

> pedigree

| ID  | Sire | Dam |
|:---:|:----:|:---:|
|  1  |  0   |  0  |
|  3  |  1   |  2  |
|  4  |  0   |  0  |
|  5  |  0   |  0  |
|  6  |  4   |  5  |
|  7  |  2   |  6  |
|  8  |  3   |  6  |

### Genotype data

The genotype should be a `data.frame` or `matrix` structure, which will be used to construct the genomic relationship matrix. The individual IDs should be row names, followed by the genotype encoding for each marker in each column. The individual IDs must match those specified in the pedigree.

In genotype encoding, `0` and `2` represent the two homozygotes, while `1` represent heterozygotes. `-9`/`NA` represents missing values.

> genotype

|     | snp1 | snp2 | snp3 | snp4 | snp5 | ... | snpN |
|:---:|:----:|:----:|:----:|:----:|:----:|:---:|:----:|
|  1  |  1   |  1   |  2   |  1   |  2   |  …  |  0   |
|  2  |  1   |  1   |  0   |  1   |  0   |  …  |  2   |
|  3  |  1   |  2   |  2   |  1   |  2   |  …  |  0   |
|  4  |  1   |  1   |  2   |  1   |  2   |  …  |  0   |
|  5  |  0   |  0   |  0   |  0   |  0   |  …  |  0   |
| ... | ...  | ...  | ...  | ...  | ...  | ... | ...  |

### Intermediate omics data

The intermediate omics should be a `data.frame` or `matrix` structure, which will be used to construct the omics similarity matrix, with $O = \frac{MM'}{N}$, where `M` is the omics feature measures (individual as rows, measures as columns), `N` represents feature numbers.

The individual IDs should be row names, followed by the omics feature measures for each feature in each column. The individual IDs must match those specified in the pedigree. Missing value should be encoding with `NA`.

> omics

|     | gene1 | gene2 | gene3 | gene4 | gene5 | ... | geneN |
|:---:|:-----:|:-----:|:-----:|:-----:|:-----:|:---:|:-----:|
|  1  | 13.86 | 15.06 | 10.89 | 14.93 | 12.86 |  …  | 14.33 |
|  2  | 0.05  | 0.06  | 0.17  | 0.11  | 0.04  |  …  | 0.08  |
|  3  | 8.36  | 8.11  | 9.49  | 10.28 | 7.322 |  …  | 7.77  |
|  4  | 0.59  | 0.43  | 1.13  | 0.51  | 0.49  |  …  | 0.48  |
|  5  | 0.24  | 0.25  | 2.45  | 0.37  | 0.14  |  …  | 0.18  |
|  6  | 11.87 | 12.42 | 9.79  | 10.78 | 8.97  | ... | 10.48 |
| ... |  ...  |  ...  |  ...  |  ...  |  ...  | ... |  ...  |

## USAGE

### Basic

To run `FSBLUP`, several basic data should provide: the phenotype data (IDs in first column) and three matrices (IDs in row and col names) or pedigree (three columns: ID\|Sire\|Dam), genomic (IDs in row names), and omics data (IDs in row names). `FSBLUP` will try to adjust same order for data before analysis.

0.  Read example data

``` r
library(tidyverse)
library(FSBLUP)
phenotype = read.csv("../test/test_pheno.csv")
pedi = read.csv("../test/test_pedi.csv")
snp = read.csv("../test/test_snp.csv") %>% column_to_rownames(colnames(.)[1]) %>% as.matrix()
omic = read.csv("../test/test_omic.csv") %>% column_to_rownames(colnames(.)[1]) %>% as.matrix()
A_mat = read.csv("../test/test_matA.csv")  %>% column_to_rownames(colnames(.)[1]) %>% as.matrix()
G_mat = read.csv("../test/test_matG.csv") %>% column_to_rownames(colnames(.)[1]) %>% as.matrix()
O_mat = read.csv("../test/test_matO.csv") %>% column_to_rownames(colnames(.)[1]) %>% as.matrix()
```

1.  Run with supplied matrices

``` r
fs = FSBLUP(phe = phenotype, trait_col = 2, M1 = A_mat, M2 = G_mat, M3 = O_mat)
# phe: phenotype
# trait_col, the column number of predicted phenotype
# M1/M2/M3: basic kinship matrix calculated by pedigree/snp/omics data
```

2.  Run with origin data

``` r
fs = FSBLUP(phe = phenotype, trait_col = 2, pedi = pedi, snp = snp, omic = omic, FPKM.qc = TRUE)
# phe: phenotype
# trait_col, the column number of predicted phenotype
# pedi: pedigree data (three columns: ID\Sire\Dam)
# snp: snp data (code as 0/1/2)
# omic: intermediate omics data
# FPKM.qc: if omic data is transcriptomics FPKM format data, will delete transcript <0.1 in more than 95% individuals.
```

3.  Change the sample number and validation number for cross_validation:

``` r
fs = FSBLUP(phe = phenotype, trait_col = 2, M1 = A_mat, M2 = G_mat, M3 = O_mat, po.crv.num = 5, po.crv.rep.num = 2)
# phe: phenotype
# trait_col, the column number of predicted phenotype
# M1/M2/M3: basic kinship matrix calculated by pedigree/snp/omics data
# po.crv.num: fold number of the cross validation in parameter optimization, default 5
# po.crv.rep.num: repeat number of the cross validation in parameter optimization default 2
```

4.  Change the classification variable for time order validation:

``` r
fs = FSBLUP(phe = phenotype, trait_col = 2, M1 = A_mat, M2 = G_mat, M3 = O_mat, po.year = 2024, po.ngen = 5)
# phe: phenotype
# trait_col, the column number of predicted phenotype
# M1/M2/M3: basic kinship matrix calculated by pedigree/snp/omics data
# po.year: the year for seprate training and validation individuals in parameter optimization, require "year" column in phenotype data
# po.ngen: the ngen for seprate training and validation individuals in parameter optimization, require "ngen" column in phenotype data
# Note should only provide one of parameters (po.year, po.ngen)
```

5.  Change the start points number of grid search procedure

``` r
fs = FSBLUP(phe = phenotype, trait_col = 2, M1 = A_mat, M2 = G_mat, M3 = O_mat, po.gs.point_num = 25)
# phe: phenotype
# trait_col, the column number of predicted phenotype
# M1/M2/M3: basic kinship matrix calculated by pedigree/snp/omics data
# po.gs.point_num: start points number, \sqrt(po.gs.point_num) is the iteration step of parameters in two direction, default 25
```

6.  Change the threshold between two iterations, and the maximum iteration number of adaptive bisection algorithm

``` r
fs = FSBLUP(phe = phenotype, trait_col = 2, M1 = A_mat, M2 = G_mat, M3 = O_mat, po.bi.max_iter = 10, po.bi.threshold = 1e-4)
# phe: phenotype
# trait_col, the column number of predicted phenotype
# M1/M2/M3: basic kinship matrix calculated by pedigree/snp/omics data
# po.bi.max_iter: the max number of iteration for bisection algorithm, default 10
# po.bi.threshold: the threshold between two iterations, stop bisection if reaches, default 1e-4
```

### Advanced

7.  Only to optimize a fusion matrix

``` r
fs = FSBLUP(phe = phenotype, trait_col = 2, M1 = A_mat, M2 = G_mat, M3 = O_mat, return.matrix = TRUE)
# phe: phenotype
# trait_col, the column number of predicted phenotype
# M1/M2/M3: basic kinship matrix calculated by pedigree/snp/omics data
# return.matrix: return the fusion matrix instead of predicted GEBVs
```

8.  Change the indicator in parameter optimization procedure

``` r
fs = FSBLUP(phe = phenotype, trait_col = 2, M1 = A_mat, M2 = G_mat, M3 = O_mat, stas.type = "rmse")
# phe: phenotype
# trait_col, the column number of predicted phenotype
# M1/M2/M3: basic kinship matrix calculated by pedigree/snp/omics data
# stas.type: indicator when optimize the fusion matrix, e.g. "cor", "auc", "rmse", "other", default as cor
```

9.  Change the column when calculating the indicator in parameter optimization procedure

``` r
fs = FSBLUP(phe = phenotype, trait_col = 2, M1 = A_mat, M2 = G_mat, M3 = O_mat, stas.phe.col = 3)
# phe: phenotype
# trait_col, the column number of predicted phenotype
# M1/M2/M3: basic kinship matrix calculated by pedigree/snp/omics data
# stas.phe.col: when calculating the indicator, use the specified column instead of the analysis column, e.g. cor(GEBV, stas.phe.col)
```

10.Change the indicator calculation with the custom function

``` r
fs = FSBLUP(phe = phenotype, trait_col = 2, M1 = A_mat, M2 = G_mat, M3 = O_mat, stas.type = "other", stas.fn = fn)
# phe: phenotype
# trait_col, the column number of predicted phenotype
# M1/M2/M3: basic kinship matrix calculated by pedigree/snp/omics data
# stas.type: indicator when optimize the fusion matrix, e.g. "cor", "auc", "rmse", "other", default as cor
# stas.fn: custom function, valid only when "type = other". 
# Note that in custom function, only supplise two parameters (X1 and X2), X1 represents the predicted GEBVs, X2 represents the target column in phentoype 
# e.g. fn = function(X1, X2) {val = cor(X1, X2); return(val)}
```

## OUTPUT

`FSBLUP` default return a list including: genetic variance (Vu), residual variance (Ve), coefficients of all fixed effects (beta), predicted GEBVs (u), maximized log-likelihood (LL)

``` r
> str(fs)
List of 5
 $ Vu  : num 3.19
 $ Ve  : num 213
 $ beta: num [1(1d)] 1.11
 $ u   : num [1:198(1d)] -0.928 -0.51 0.421 0.378 0.449 ...
  ..- attr(*, "dimnames")=List of 1
  .. ..$ : chr [1:198] "754" "986" "993" "994" ...
 $ LL  : num -809
```

Once use `return.matrix = TRUE` parameters, the returns will change to:

``` r
> fs[1:10, 1:10]
      11211 11410 11413 11553 11590 11625 11643 11649 11685 11715
11211  0.90  0.17  0.17  0.02  0.17  0.19  0.04  0.17  0.07 -0.26
11410  0.17  1.59  0.99  0.17 -0.10 -0.21 -0.17  0.73  0.69 -1.00
11413  0.17  0.99  1.05  0.11 -0.08  0.00 -0.08  0.40  0.31 -0.53
11553  0.02  0.17  0.11  0.60 -0.01 -0.13 -0.07  0.02 -0.04 -0.10
11590  0.17 -0.10 -0.08 -0.01  0.94  0.08 -0.10  0.18  0.17 -0.03
11625  0.19 -0.21  0.00 -0.13  0.08  0.96 -0.04 -0.15 -0.25  0.37
11643  0.04 -0.17 -0.08 -0.07 -0.10 -0.04  0.59 -0.09 -0.18  0.08
11649  0.17  0.73  0.40  0.02  0.18 -0.15 -0.09  1.16  1.03 -0.80
11685  0.07  0.69  0.31 -0.04  0.17 -0.25 -0.18  1.03  1.48 -0.96
11715 -0.26 -1.00 -0.53 -0.10 -0.03  0.37  0.08 -0.80 -0.96  1.98
```
