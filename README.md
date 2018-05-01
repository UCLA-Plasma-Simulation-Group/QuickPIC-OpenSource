# ![QUICKPIC](http://exodus.physics.ucla.edu/~uclapic/repo_images/quickpic.png)
QuickPIC is a 3D parallel (MPI & OpenMP Hybrid) Quasi-Static PIC code, which is developed based on the framework UPIC. This is the UCLA Plasma Simulation Group's official open-source repository for QuickPIC.

# Upon cloning the repository

If you clone this repository, we ask that you __please contact__ Weiming An (anweiming@ucla.edu). The development of QuickPIC relies on grant funding for which we are required to report code usage, and we would greatly appreciate being able to maintain accurate user-number tallies.

Please also feel free to send an email to this address if you would like assistance using the code, if you have any questions related to the code base, or if you have suggestions for improvements. We are in the process of establishing a user forum for discussions related to the use and development of QuickPIC and will be more than happy to include you in such a forum if you desire.

# Compile QuickPIC

The makefile is set to use gcc and gfortran with MPI. HDF5-Parallel is also required
for compiling. 

To compile the programs, execute:
```
make
```
The program name is qpic.e

The command to execute a program with both MPI and OpenMP varies from
one system to another.  One possible command is:
```
mpirun -np nproc ./qpic.e
```
where nproc is the number of processors to be used. Note that nproc should be 2 at least.

By default, OpenMP will use the maximum number of processors it can find
on the MPI node.  If the user wants to control the number of threads, the
environment variable 'OMP_NUM_THREADS' may need to be set to the maximum
number of threads per node expected.

## Compile on Debian (Ubuntu)

To compile QuickPIC on Debian, execute:
```
make SYS_FX=DEBIAN
```
This requires that libhdf5-openmpi-dev is installed as it uses the h5pfc compiler wrapper.

# Generate documentation

QuickPIC uses [FORD](https://github.com/cmacmackin/ford) to generate documentation. Before installing FORD, make sure Python and relevant libraries are correctly installed. The simplest way to install FORD is using [pip](https://pip.pypa.io/en/latest/) with this command
```
sudo pip install ford
```
which will install FORD and all its [dependencies](https://github.com/cmacmackin/ford/wiki/Dependencies) automatically. If you wish to install the most recent version (in case that FORD installed through pip runs abnormally in some systems), you can also clone the github repository:
```
git clone https://github.com/cmacmackin/ford.git
```
or download the [install package](https://github.com/cmacmackin/ford/archive/master.zip). Go into the FORD folder and run the command below to install manually:
```
python install setup.py
```

The project file (FORD_DOC.md) includes various [options](https://github.com/cmacmackin/ford/wiki/Project-File-Options) and information of the QuickPIC documentation in the meta-data. Be aware that the option __include__, which defines the directories for file searching of Fortran's intrinsic `include` statement, should be properly configured before generating the documentation. Keep in mind that you must make sure the QuickPIC source code compiles correctly because FORD is not made to check any syntax error. Run the command below to generate the documentation:
```
ford FORD_DOC.md
```