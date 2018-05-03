src_dir: ./source
output_dir: ./doc
project: QuickPIC
project_github: https://github.com/UCLA-Plasma-Simulation-Group/QuickPIC-OpenSource
title: QuickPIC
summary: A full 3D quasi-static particle-in-cell code
author: Weiming An
author_description: 
github: https://github.com/UCLA-Plasma-Simulation-Group/QuickPIC-OpenSource
email: anweiming@ucla.edu
fpp_extensions: fpp
predocmark: >
media_dir: ./media
docmark_alt: #
predocmark_alt: <
display: public
         protected
         private
source: false
graph: true
search: true
macro: TEST
       LOGIC=.true.
extra_mods: json_module: http://jacobwilliams.github.io/json-fortran/
            futility: http://cmacmackin.github.io
license: by-nc
extra_filetypes: sh #
include: /usr/include/mpi

QuickPIC is a 3D parallel (MPI & OpenMP Hybrid) Quasi-Static PIC code, which is developed based on the framework UPIC. This is the UCLA Plasma Simulation Group's official open-source repository for QuickPIC.

# Upon cloning the repository

If you clone this repository, we ask that you please contact Weiming An (anweiming@ucla.edu). The development of QuickPIC relies on grant funding for which we are required to report code usage, and we would greatly appreciate being able to maintain accurate user-number tallies.

Please also feel free to send an email to this address if you would like assistance using the code, if you have any questions related to the code base, or if you have suggestions for improvements. We are in the process of establishing a user forum for discussions related to the use and development of QuickPIC and will be more than happy to include you in such a forum if you desire.

# Compile QuickPIC

The makefile is set to use gcc and gfortran with MPI. HDF5-Parallel is also required for compiling.

To compile the programs, execute:
```
make
```
The program name is qpic.e

The command to execute a program with both MPI and OpenMP varies from one system to another. One possible command is:
```
mpirun -np nproc ./qpic.e
```
where nproc is the number of processors to be used. Note that nproc should be 2 at least.

By default, OpenMP will use the maximum number of processors it can find on the MPI node. If the user wants to control the number of threads, the environment variable 'OMP_NUM_THREADS' may need to be set to the maximum number of threads per node expected.

# Compile on Debian (Ubuntu)

To compile QuickPIC on Debian, execute:
```
make SYS_FX=DEBIAN
```
This requires that libhdf5-openmpi-dev is installed as it uses the h5pfc compiler wrapper.

<!-- @Note -->
<!-- You can include any notes (or bugs, warnings, or todos) like so. -->

<!-- @Bug
No bug found in current version.
@endbug -->

