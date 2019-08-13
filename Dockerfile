from bioconductor/devel_core2
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

ENV BIO_VERSION devel

RUN echo 'local({\nr <- getOption("repos")\nr["BioC"] <- "https://bioconductor.org/packages/${BIO_VERSION}/bioc"\nr["CRAN"] <- "http://cran.rstudio.com/"\noptions(repos = r)\nr["BioCann"] <- "https://bioconductor.org/packages/${BIO_VERSION}/data/annotation"\nr["BioCexp"] <- "https://bioconductor.org/packages/${BIO_VERSION}/data/experiment"\n})' > ~/.Rprofile \
    && Rscript -e "install.packages(c('devtools', 'BiocManager', 'remotes'), dependencies=TRUE, ask = FALSE)" \
    # NOTE: these seem to be required for garnett to succeed in docker
    && Rscript -e "BiocManager.install(c('org.Hs.eg.db', 'org.Mm.eg.db', 'HSMMSingleCell'), dependencies=TRUE, ask = FALSE)" \
    && Rscript -e "devtools::install_github(repo = 'bimberlabinternal/OOSAP', ref = 'SingleR', dependencies = T, upgrade = 'always', ask=FALSE)" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds