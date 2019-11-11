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

# Let this run for the purpose of installing/caching dependencies
RUN Rscript -e "install.packages(c('devtools', 'BiocManager', 'remotes'), dependencies=TRUE, ask = FALSE)" \
    && echo -e "local({\noptions(repos = BiocManager::repositories())\n})\n" >> ~/.Rprofile.site \
    # NOTE: these seem to be required for garnett to succeed in docker. DESeq2/genefilter added for Seurat
    && Rscript -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db', 'HSMMSingleCell', 'monocle', 'DelayedMatrixStats', 'DESeq2', 'genefilter'), dependencies=TRUE, ask = FALSE)" \
    && Rscript -e "devtools::install_github(repo = 'bimberlabinternal/OOSAP', ref = 'Dev', dependencies = T, upgrade = 'always')" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# This should not be cached if the files change
ADD . /OOSAP

RUN cd /OOSAP \
    && R CMD build . \
    && R CMD INSTALL --build *.tar.gz \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds