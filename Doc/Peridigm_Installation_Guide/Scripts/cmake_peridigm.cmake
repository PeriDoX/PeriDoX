rm -f CMakeCache.txt
rm -rf CMakeFiles/

cmake \
-D CMAKE_BUILD_TYPE:STRING=Release \
-D Trilinos_DIR:PATH=/usr/local/bin/trilinos-12.4.2/lib/cmake/Trilinos/ \
-D CMAKE_C_COMPILER:STRING=/usr/local/lib/openmpi-1.10.2/bin/mpicc \
-D CMAKE_CXX_COMPILER:STRING=/usr/local/lib/openmpi-1.10.2/bin/mpicxx \
-D BOOST_ROOT=/usr/local/lib/boost-1.55.0/ \
-D CMAKE_CXX_FLAGS:STRING="-O2 -Wall -std=c++11 -pedantic -Wno-long-long -ftrapv -Wno-deprecated" \
/usr/local/src/Peridigm-1.4.1-Source