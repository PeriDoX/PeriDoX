#!/bin/bash

#######################################
# Header                              #
#######################################
#
# Install Peridigm from github repository master
#
# Requirements:
#            Allocate exclusive node before using this script:
#            salloc --exclusive
# 
# Revisions: 2017-12-22 Martin Raedel <martin.raedel@dlr.de>
#                       Initial draft
#               
# Contact:   Martin Raedel,  martin.raedel@dlr.de
#            DLR Composite Structures and Adaptive Systems
#          
#                                 __/|__
#                                /_/_/_/  
#            www.dlr.de/fa/en      |/ DLR
# 
#######################################
# Content                             #
#######################################

#--------------------------------------
# Variables
#--------------------------------------

# The directory where all the magic happens - should not exist in advance
basedir=$HOME'/Documents/Peridigm/20171222EC_/'

# Path from where to clone
githubclonepath='https://github.com/peridigm/peridigm.git'

# Internal directory names
builddir='build'
srcdir='src'

# File names
file_cmake='cmake_peridigm.cmake'
file_cmake_log='cmake_peridigm.log'
file_make_log='make.log'
file_make_test_log='make_test.log'
file_make_install_log='make_install.log'

file_peridigm_bin='Peridigm'

# Number of CPUs for make
make_cpus = 8

#--------------------------------------
# Script
#--------------------------------------

#------------------
# Load build environment
#------------------

echo 'Load build environment'
. /cluster/software/slurm/etc/env.d/mpibuild.sh
. /cluster/software/slurm/etc/env.d/peridigm.sh

#------------------
# Folder structure
#------------------

echo 'Create directory structure'
if [ -d ${basedir} ]; then  # Control will enter here if $DIRECTORY doesn't exist.
  echo 'Directory '${basedir}' already exists. Exit.'
  exit 0
fi

mkdir ${basedir}
cd ${basedir}
mkdir ${builddir}
mkdir ${srcdir}
cd ${srcdir}

#------------------
# Clone from GitHub
#------------------

echo 'Clone from GitHub'
git clone ${githubclonepath}
cd ../${builddir}

#------------------
# Create cmake file
#------------------

echo 'Create and execute cmake file'

# Fill file
echo 'rm -f CMakeCache.txt' >> ${file_cmake}
echo 'rm -rf CMakeFiles/' >> ${file_cmake}
echo '' >> ${file_cmake}
echo 'cmake \' >> ${file_cmake}
echo '-D CMAKE_BUILD_TYPE:STRING=Release \' >> ${file_cmake}
echo '-D CMAKE_INSTALL_PREFIX='${basedir}${builddir}' \' >> ${file_cmake}
echo '-D CMAKE_CXX_FLAGS:STRING="-O2 -Wall -std=c++11 -pedantic -Wno-long-long -ftrapv -Wno-deprecated" \' >> ${file_cmake}
echo ${basedir}${srcdir}'/peridigm/' >> ${file_cmake}

# Change permissions to make cmake-file executable
chmod u+x ${file_cmake}

# Execute cmake-file
./${file_cmake} > ${file_cmake_log} 2>&1

#------------------
# Make
#------------------

echo 'make'
make -j ${make_cpus} > ${file_make_log} 2>&1

echo 'make test'
make test > ${file_make_test_log} 2>&1 
if grep -q Error ${file_make_test_log}; then
  echo '  make test contains failed tests'
else
  echo '  all tests passed'
fi

echo 'make install'
make install > ${file_make_install_log} 2>&1

#------------------
# Comment
#------------------

cd bin
echo 'Finished - Peridigm executable path':
readlink -f ${file_peridigm_bin}

#--------------------------------------
# Clean
#--------------------------------------