#!/bin/bash

# See:
# http://dirk.eddelbuettel.com/blog/2017/11/27/
# https://www.jumpingrivers.com/blog/speeding-up-package-installation/
# https://github.com/travis-ci/travis-build/blob/master/lib/travis/build/script/r.rb
echo 'Running travis setup'

mkdir $HOME/.R

export R_BIOC_VERSION='3.10'

# Force python3:
python3 -m pip install --upgrade pip

python3 --version
python --version
pip3 --version
pip --version

echo -e 'CXX_STD = CXX14\n\nVER=\nCCACHE=ccache\nCC=$(CCACHE) gcc$(VER) -std=gnu99\nCXX=$(CCACHE) g++$(VER)\nC11=$(CCACHE) g++$(VER)\nC14=$(CCACHE) g++$(VER)\nFC=$(CCACHE) gfortran$(VER)\nF77=$(CCACHE) gfortran$(VER)\n' > $HOME/.R/Makevars
echo 'max_size = 5.0G\n# important for R CMD INSTALL *.tar.gz as tarballs are expanded freshly -> fresh ctime\nsloppiness = include_file_ctime\n# also important as the (temp.) directory name will differ\nhash_dir = false' > ~/.ccache/ccache.conf

# See: https://stackoverflow.com/questions/49525561/rcppeigen-package-pragma-clang-diagnostic-pop-warnings
echo -e 'PKG_CXXFLAGS=-w -DEIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS' >> $HOME/.R/Makevars

CORES=`Rscript -e "getOption('Ncpus', 1L)"`
echo "Existing Ncpus: $CORES"

CORES=`Rscript -e "parallel::detectCores()"`
echo "Detected Cores: $CORES"

CORES=2
echo "options(Ncpus = ${CORES})" >> ~/.Rprofile
CORES=`Rscript -e "getOption('Ncpus', 1L)"`
echo "Final Ncpus: $CORES"

# Verify config, add Bioconductor repos to default repos, so install.packages will pick them up.
if [ -e ~/.Rprofile ];then
    echo 'Rprofile exists:'
    cat ~/.Rprofile
fi

if [ -e ~/.Rprofile.site ];then
    echo 'Rprofile.site exists:'
    cat ~/.Rprofile.site
fi

Rscript -e "install.packages(c('devtools', 'BiocManager', 'remotes'), dependencies=TRUE, ask = FALSE)"

# See: https://stackoverflow.com/questions/26042751/cannot-install-package-xml-to-r
echo "options(repos = c(BiocManager::repositories(version = '${R_BIOC_VERSION}'), c('xml' = 'http://www.omegahat.net/R')));" >> ~/.Rprofile.site

# Log repos to ensure Bioconductor used:
echo 'R repos:'
Rscript -e "getOption('repos')"

# Bioconductor
Rscript -e "BiocManager::install(version='${R_BIOC_VERSION}')"

# This is to pre-install two packages with exceptionally large logging, to avoid travis-ci log size errors:
trap 'catch' ERR
catch() {
  cat ~/installWithoutConsole.log
  exit 1
}

echo 'Installing packages with output suppressed to reduce log size'
Rscript -e "install.packages(c('RSQLite', 'igraph', 'RcppEigen', 'glmnet', 'rJava'), dependencies=TRUE, ask = FALSE)" > ~/installWithoutConsole.log
echo 'Installing more packages with output suppressed to reduce log size'
Rscript -e "install.packages(c('BiocParallel', 'viridis', 'isoband', 'Rgraphviz'), dependencies=TRUE, ask = FALSE)" >> ~/installWithoutConsole.log
echo 'Done installing'