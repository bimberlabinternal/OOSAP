from bioconductor/release_core2
RUN apt-get update -y \
	&& apt-get install -y \
		libhdf5-dev \
		python-pip \
    && pip install wheel \
    && pip install numba==0.42.0 \
    && pip install umap-learn \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Note: SingleR depends on Seurat, but doesnt actually install it.  For the time being, manually install Seurat here; however,
# Once SingleR improves we can hopefully remove this and let the normal R installation take care of it
RUN Rscript -e "install.packages(c('devtools', 'remotes', 'Seurat'), dependencies=TRUE, repos='http://cran.rstudio.com/')" \
    && Rscript -e "BiocManager::install(c('monocle', 'AnnotationDbi', 'Biobase', 'DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db', 'loomR', 'SingleCellExperiment', 'MAST', 'DESeq2'))" \
    && Rscript -e "install.packages(c('Seurat'), dependencies=TRUE, repos='http://cran.rstudio.com/')" \
    && Rscript -e "devtools::install_github(repo = 'bimberlabinternal/OOSAP', ref = 'Dev', dependencies = T, upgrade = 'always', ask=FALSE)" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds