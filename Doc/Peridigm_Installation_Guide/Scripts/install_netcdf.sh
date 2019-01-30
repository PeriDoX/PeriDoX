# Set environment variables for MPI compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77

# Configure NetCDF
./configure --prefix=/usr/local/bin/netcdf-4.4.0/ --disable-netcdf-4 --disable-dap

# Make, test, and install NetCDF
make -j 4
make check
make install