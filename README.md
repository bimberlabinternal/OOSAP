[![Build Status](https://travis-ci.com/bimberlabinternal/OOSAP.svg?branch=Dev)](https://travis-ci.com/bimberlabinternal/OOSAP)

# OOSAP

OHSU-ONPRC Single-cell Analysis Package

The goal of this package is bringing together tools and pipelines for the analysis of single-cell RNASeq data.

### Installation

```{r }

# Install requirements.  Other dependencies will be downloaded automatically
install.packages("devtools", 'remotes')

#Latest version:
install_github("bimberlabinternal/OOSAP")

#Or a specific release :
install_github("bimberlabinternal/OOSAP", tag = "1.0")

```
Pre-packaged Docker images with all needed dependencies installed can be found on our [dockerhub repository](https://hub.docker.com/r/bimberlab/oosap): 

```

docker pull bimberlab/oosap

```
