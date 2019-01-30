# Set environment variables for MPI compilers
export CC=mpicc
export CXX=mpicxx
export FC=mpif90
export F77=mpif77

# Run the Boost bootstrap script
./bootstrap.sh

# add using mpi to project-config.jam
echo "using mpi ;" >> project-config.jam

# Compile and install Boost using the Boost's bjam build system
./b2 install --prefix=/usr/local/lib/boost-1.60.0/
