#!/bin/bash

# See:
# http://dirk.eddelbuettel.com/blog/2017/11/27/
# https://www.jumpingrivers.com/blog/speeding-up-package-installation/
echo 'Running travis setup'

mkdir $HOME/.R
Rscript -e "getOption("Ncpus", 1L)"
CORES=`Rscript -e "parallel::detectCores()"`
echo "Cores: $CORES"

echo -e 'CXX_STD = CXX14\n\nVER=\nCCACHE=ccache\nCC=$(CCACHE) gcc$(VER) -std=gnu99\nCXX=$(CCACHE) g++$(VER)\nC11=$(CCACHE) g++$(VER)\nC14=$(CCACHE) g++$(VER)\nFC=$(CCACHE) gfortran$(VER)\nF77=$(CCACHE) gfortran$(VER)\nMAKE=make -j4' > $HOME/.R/Makevars
echo "options(Ncpus = ${CORES})" >> ~/.Rprofile
echo 'max_size = 5.0G\n# important for R CMD INSTALL *.tar.gz as tarballs are expanded freshly -> fresh ctime\nsloppiness = include_file_ctime\n# also important as the (temp.) directory name will differ\nhash_dir = false' > ~/.ccache/ccache.conf