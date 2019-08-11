from bioconductor/release_core2
RUN apt-get update -y \
	&& apt-get install -y \
		libhdf5-dev \
		libpython-dev \
		python-pip \
    && pip install wheel \
    && pip install numba==0.42.0 \
    && pip install umap-learn \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('devtools', 'remotes'), dependencies=TRUE, repos='http://cran.rstudio.com/')" \
    # NOTE: manually installing might be required if remotes goes out of order.  monocle is needed before garnett
    # && Rscript -e "BiocManager::install(c('monocle', 'AnnotationDbi', 'Biobase', 'DelayedArray', 'DelayedMatrixStats', 'org.Hs.eg.db', 'org.Mm.eg.db', 'loomR', 'SingleCellExperiment', 'MAST', 'DESeq2'))" \
    && Rscript -e "devtools::install_github(repo = 'bimberlabinternal/OOSAP', ref = 'Dev', dependencies = T, upgrade = 'always', ask=FALSE)" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds