from bioconductor/release_core2
RUN apt-get update -y \
	&& apt-get install -y \
		libhdf5-dev \
		python-pip \
    && pip install numba==0.42.0
    && pip install umap-learn
    && apt-get remove -y python-pip
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

RUN Rscript -e "install.packages(c('devtools', 'remotes') ,dependencies=TRUE, repos='http://cran.rstudio.com/')" &&
    Rscript -e "devtools::install_github(repo = 'bimberlabinternal/OOSAP', ref = 'Dev', dependencies = T, upgrade = 'always', ask=FALSE)" &&
    rm -rf /tmp/downloaded_packages/ /tmp/*.rds