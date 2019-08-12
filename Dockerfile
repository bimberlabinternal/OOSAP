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

RUN echo 'local({\nr <- getOption("repos")\nr["BioC"] <- "https://bioconductor.org/packages/release/bioc"\nr["CRAN"] <- "http://cran.rstudio.com/"\noptions(repos = r)\nr["BioCann"] <- "https://bioconductor.org/packages/release/data/annotation"\n})' > ~/.Rprofile \
    && Rscript -e "install.packages(c('devtools', 'BiocManager', 'remotes'), dependencies=TRUE)" \
    && Rscript -e "devtools::install_github(repo = 'bimberlabinternal/OOSAP', ref = 'Dev', dependencies = T, upgrade = 'always', ask=FALSE)" \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds