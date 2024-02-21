# BREW3R.r
The BREW3R.r package has been written to be part of the BREW3R workflow. Today, the package contains a single function which enable to extend three prime of gene annotations using another gene annotation as template. This is very helpful when you are using a technique that only sequence three-prime end of genes like 10X scRNA-seq or BRB-seq.

## Installation

To install from Bioconductor use:

```{r installation bioconductor, eval=FALSE}
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("BREW3R.r")
```

To install from github use:

```{r installation github, eval=FALSE}
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("lldelisle/BREW3R.r")
```

## Issues
If you have issues, use the Issues in github or send an email to lucille.delisle\@epfl.ch
