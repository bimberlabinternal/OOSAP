![R Build and Checks](https://github.com/bimberlabinternal/OOSAP/workflows/R%20Build%20and%20Checks/badge.svg)
[![DOI](https://zenodo.org/badge/187915110.svg)](https://zenodo.org/badge/latestdoi/187915110)

# OOSAP

OHSU-ONPRC Single-cell Analysis Package

The goal of this package is to provide a high-level wrapper around tools and pipelines for the analysis of single-cell RNASeq data.

## Table of Contents
* [Installation](#installation)
* [Development Guidelines](#developers)

### <a name="installation">Installation</a>

```{r }

# Install requirements.  Other dependencies will be downloaded automatically
install.packages(c('devtools', 'BiocManager', 'remotes'), dependencies=TRUE, ask = FALSE)

# Make sure to update your Rprofile i.e., ~/.Rprofile.site
# local({options(repos = BiocManager::repositories())})

# Get some of the core packages
BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db', 'HSMMSingleCell', 'monocle', 'DelayedMatrixStats', 'DESeq2', 'genefilter'), dependencies=TRUE, ask = FALSE)

#Latest version:
devtools::install_github(repo = 'bimberlabinternal/OOSAP', ref = 'Dev', dependencies = T, upgrade = 'always')


```
    
    
Pre-packaged Docker images with all needed dependencies installed can be found on our [GitHub Packages repository](https://github.com/orgs/bimberlabinternal/packages/container/package/oosap): 

```

docker pull ghcr.io/bimberlabinternal/oosap:latest

```

### <a name="developers">Development Guidelines</a>

* New development should occur on a branch, and go through a Pull Request before merging into Dev.  [See here for information on the pull request workflow](https://guides.github.com/introduction/flow/).  Ideally PRs would be reviewed by another person.  For the PR, please review the set of changed files carefully to make sure you are only merging the changes you intend.   

* New exported functions should have [Roxygen2 documentation](https://kbroman.org/pkg_primer/pages/docs.html).

* As part of each PR, you should run 'devtools::document()' to update documentation and include these changes with your commits.

* It is a good idea to run 'R CMD check' locally to make sure your changes will pass.  [See here for more information](http://r-pkgs.had.co.nz/check.html)

* Code should only be merged after the Travis CI build and tests pass.  The Dev branch should always be stable.

* New features should ideally have at least a basic test (see [R testthat](http://r-pkgs.had.co.nz/tests.html)).  There is existing test data in ./tests/testdata.  This can be expanded, but please be conscious about file size and try to reuse data across tests if appropriate.

  
