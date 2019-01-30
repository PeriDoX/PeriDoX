# Set environment variables for MPI compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77

# Configure HDF5
./configure --prefix=/usr/local/bin/hdf5-1.8.16/ --enable-parallel

# Make and install HDF5
make -j 4
make install
