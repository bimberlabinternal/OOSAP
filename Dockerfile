from bioconductor/bioconductor_docker:RELEASE_3_10

# NOTE: if anything breaks the dockerhub build cache, you will probably need to build locally and push to dockerhub.
# After the cache is in place, builds from github commits should be fast.
RUN apt-get update -y \
	&& apt-get upgrade -y \
	&& apt-get install -y \
		libhdf5-dev \
		libpython3-dev \
		python3-pip \
		openjdk-11-jdk \
	&& pip3 install umap-learn \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*

# Let this run for the purpose of installing/caching dependencies
RUN Rscript -e "install.packages(c('devtools', 'BiocManager', 'remotes'), dependencies=TRUE, ask = FALSE)" \
	&& echo "local({\noptions(repos = BiocManager::repositories())\n})\n" >> ~/.Rprofile \
	# See: https://stackoverflow.com/questions/26042751/cannot-install-package-xml-to-r
	&& Rscript -e "install.packages('XML', repos = 'http://www.omegahat.net/R');" \
	# NOTE: these seem to be required for garnett to succeed in docker. DESeq2/genefilter added for Seurat
	&& Rscript -e "BiocManager::install(c('org.Hs.eg.db', 'org.Mm.eg.db', 'HSMMSingleCell', 'monocle', 'DelayedMatrixStats', 'DESeq2', 'genefilter'), dependencies=TRUE, ask = FALSE)" \
    && Rscript -e "devtools::install_github(repo = 'bimberlabinternal/OOSAP', ref = 'Dev', dependencies = T, upgrade = 'always')" \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# This should not be cached if the files change
ADD . /OOSAP

RUN cd /OOSAP \
	&& R CMD build . \
	&& Rscript -e "print(getOption('repos'))" \
	&& Rscript -e "BiocManager::install(ask = F);" \
	&& Rscript -e "devtools::install_deps(pkg = '.', dependencies = TRUE, upgrade = 'always');" \
	&& R CMD INSTALL --build *.tar.gz \
	&& rm -Rf /tmp/downloaded_packages/ /tmp/*.rds