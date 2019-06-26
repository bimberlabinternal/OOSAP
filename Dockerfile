from bioconductor/release_core2
RUN apt-get update -y \
	&& apt-get install -y \
		libhdf5-dev \
		python-pip \
	&& apt-get clean \
	&& rm -rf /var/lib/apt/lists/*
RUN pip install numba==0.42.0 && pip install umap-learn && apt-get remove -y python-pip

ADD scripts/install.R /
RUN Rscript install.R && rm -rf /tmp/downloaded_packages/ /tmp/*.rds