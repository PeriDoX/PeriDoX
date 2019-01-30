# Set environment variables for MPI compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77

# Run the Boost bootstrap script
./bootstrap.sh

# add using mpi to user-config.jam
echo "using mpi ;" >> tools/build/v2/user-config.jam
cp tools/build/v2/user-config.jam ~/

# Compile and install Boost using the Boost's bjam build system
./b2 install --prefix=/usr/local/lib/boost-1.55.0/
