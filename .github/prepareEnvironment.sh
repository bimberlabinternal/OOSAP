#!/bin/bash

# See:
# http://dirk.eddelbuettel.com/blog/2017/11/27/
# https://www.jumpingrivers.com/blog/speeding-up-package-installation/

mkdir $HOME/.R

python3 --version

# ccache / R:
echo -e 'CXX_STD = CXX14\nVER=\nCCACHE=ccache\nCC=$(CCACHE) gcc$(VER) -std=gnu99\nCXX=$(CCACHE) g++$(VER)\nC11=$(CCACHE) g++$(VER)\nC14=$(CCACHE) g++$(VER)\nFC=$(CCACHE) gfortran$(VER)\nF77=$(CCACHE) gfortran$(VER)' > $HOME/.R/Makevars

if [ ! -e ~/.ccache/ ];then
    mkdir -p ~/.ccache/
fi
echo 'max_size = 5.0G\n# important for R CMD INSTALL *.tar.gz as tarballs are expanded freshly -> fresh ctime\nsloppiness = include_file_ctime\n# also important as the (temp.) directory name will differ\nhash_dir = false' > ~/.ccache/ccache.conf

# See: https://stackoverflow.com/questions/49525561/rcppeigen-package-pragma-clang-diagnostic-pop-warnings
addFlags() {
    VAR=$1

    VALUE=`R CMD config ${VAR}`
    echo $VAR": "$VALUE
    if [[ ${VALUE} != *"ERROR"* ]];then
        echo -e $VAR"="${VALUE}" -w -pipe -pedantic -DEIGEN_PERMANENTLY_DISABLE_STUPID_WARNINGS" >> $HOME/.R/Makevars
        R CMD config ${VAR}
    fi
}

addFlags 'CXXFLAGS'
addFlags 'CXX1XFLAGS'
addFlags 'CXX11FLAGS'
addFlags 'CXX14FLAGS'
addFlags 'CXX17FLAGS'

echo 'Makevars file:'
cat $HOME/.R/Makevars

CORES=`Rscript -e "getOption('Ncpus', 1L)"`
echo "Existing Ncpus: $CORES"

CORES=`Rscript -e "parallel::detectCores()"`
echo "Detected Cores: $CORES"

CORES=2
echo "options(Ncpus = ${CORES})" >> ~/.Rprofile
CORES=`Rscript -e "getOption('Ncpus', 1L)"`
echo "Final Ncpus: $CORES"