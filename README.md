# ![QUICKPIC](http://exodus.physics.ucla.edu/~uclapic/repo_images/quickpic.png)
QuickPIC is a 3D parallel (MPI & OpenMP Hybrid) Quasi-Static PIC code, which is developed based on the framework UPIC. This is the UCLA Plasma Simulation Group's official open-source repository for QuickPIC.

# Upon cloning the repository

If you clone this repository, we ask that you __please contact__ Weiming An (anweiming@ucla.edu). The development of QuickPIC relies on grant funding for which we are required to report code usage, and we would greatly appreciate being able to maintain accurate user-number tallies.

Please feel free to join in our email list [quickpic-users-developers](https://groups.google.com/a/lists.ucla.edu/forum/#!forum/quickpic-users-developers). We also have created a chat room “QuickPIC-User-Developer-Community” on [Gitter](https://gitter.im/home).

# Compile QuickPIC

The makefile is set to use gcc and gfortran with MPI. HDF5-Parallel is also required
for compiling. 

To compile the programs, execute:

make

The program name is qpic.e

The command to execute a program with both MPI and OpenMP varies from
one system to another.  One possible command is:

mpirun -np nproc ./qpic.e

where nproc is the number of processors to be used. Note that nproc should be 2 at least.

By default, OpenMP will use the maximum number of processors it can find
on the MPI node.  If the user wants to control the number of threads, the
environment variable 'OMP_NUM_THREADS' may need to be set to the maximum
number of threads per node expected.

## Compile on Debian (Ubuntu)

To compile QuickPIC on Debian, execute:

make SYS_FX=DEBIAN

This requires that libhdf5-openmpi-dev is installed as it uses the h5pfc compiler wrapper.

