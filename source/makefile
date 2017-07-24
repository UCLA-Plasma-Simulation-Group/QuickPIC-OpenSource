SYS_SFX=GF_OPENMPI

# precision
# SGL = single precision, float(4)
# DBL = double precision, float(8)
PRECISION = DBL

INC_FILE_NAME = make.${SYS_SFX}
include ./${INC_FILE_NAME} 

# name of the compiler 
FC = ${FC_${SYS_SFX}}
CC = ${CC_${SYS_SFX}}
LINKER = ${LINKER_${SYS_SFX}}

# compiler options

OPTS_C = ${OPTS_BASE} ${OPTS_${PRECISION}_${SYS_SFX}} 

SPECIAL_OPTS= ${SPECIAL_OPTS_${SYS_SFX}}

FORMAT_FREE = ${FORMAT_FREE_${SYS_SFX}}
FORMAT_FIXED = ${FORMAT_FIXED_${SYS_SFX}}

# Objects list

OBJS_BASE = \
parallel_class.o parallel_pipe_class.o perrors_class.o dtimer.o spect2d_class.o \
spect3d_class.o \
hdf5io_class.o \
ufield2d_class.o ufield2d_lib77.o ufield2d_lib.o \
fft2d_lib.o fft2d_lib77.o fft2d_class.o \
fpois2d_lib.o fpois2d_lib77.o fpois2d_class.o \
part2d_lib.o part2d_lib77.o part2d_class.o fdist2d_class.o \
ufield3d_class.o ufield3d_lib77.o ufield3d_lib.o \
part3d_lib.o part3d_lib77.o part3d_class.o fdist3d_class.o \
field2d_class.o field3d_class.o species2d_class.o beam3d_class.o \
input_class.o


OBJS_MAIN = main.o  

OBJS = ${OBJS_BASE} ${OBJS_MAIN}

# hdf libraries 
HDF_INCLUDE_PATH = ${HDF_INCLUDE_PATH_${SYS_SFX}}
HDF_LIBPATH = ${HDF_LIBPATH_${SYS_SFX}}
HDF_LIBS = ${HDF_LIBPATH} 

# other libs
OTHER_LIBS = ${OTHER_LIBS_${SYS_SFX}}

# memory options
MEM_OPTIONS = ${MEM_OPTIONS_${SYS_SFX}}

# options for linking
# use static libs
STATIC_LINK = ${STATIC_LINK_${SYS_SFX}}
# use shared libs
SHARED_LINK = ${SHARED_LINK_${SYS_SFX}}

OPTS = ${HDF_INCLUDE_PATH} ${OPTS_C}

# Linkage rule
main :: ${OBJS} 
	${LINKER} ${MEM_OPTIONS} ${STATIC_LINK} ${SHARED_LINK} ${OBJS} ${FFTW_LIBS} ${HDF_LIBS} ${OTHER_LIBS} -o qpic.e

clean ::
	rm *.o *.mod; rm qpic.e

# Compilation rules
input_class.o : input_class.f spect3d_class.o\
                  perrors_class.o parallel_pipe_class.o parallel_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c input_class.f -o input_class.o

beam3d_class.o : beam3d_class.f fdist3d_class.o field3d_class.o spect3d_class.o \
                  perrors_class.o parallel_pipe_class.o part3d_class.o hdf5io_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c beam3d_class.f -o beam3d_class.o


species2d_class.o : species2d_class.f spect2d_class.o perrors_class.o \
                  parallel_pipe_class.o field2d_class.o hdf5io_class.o \
                  fdist2d_class.o part2d_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c species2d_class.f -o species2d_class.o


field3d_class.o : field3d_class.f spect3d_class.o perrors_class.o \
                  parallel_pipe_class.o ufield3d_class.o hdf5io_class.o 
	${FC} ${OPTS} ${FORMAT_FREE} -c field3d_class.f -o field3d_class.o

field2d_class.o : field2d_class.f spect2d_class.o perrors_class.o \
                  parallel_pipe_class.o ufield2d_class.o fft2d_class.o \
                  fpois2d_class.o hdf5io_class.o field3d_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c field2d_class.f -o field2d_class.o

part3d_class.o : part3d_class.f fdist3d_class.o ufield3d_class.o spect3d_class.o \
                  perrors_class.o parallel_pipe_class.o part3d_lib.o hdf5io_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c part3d_class.f -o part3d_class.o

part3d_lib77.o : part3d_lib77.f
	${FC} ${OPTS} ${FORMAT_FIXED} -c part3d_lib77.f -o part3d_lib77.o

part3d_lib.o : part3d_lib.f
	${FC} ${OPTS} ${FORMAT_FREE} -c part3d_lib.f -o part3d_lib.o

fdist3d_class.o : fdist3d_class.f spect3d_class.o perrors_class.o \
                   parallel_pipe_class.o 
	${FC} ${OPTS} ${FORMAT_FREE} -c fdist3d_class.f -o fdist3d_class.o

ufield3d_lib77.o : ufield3d_lib77.f
	${FC} ${OPTS} ${FORMAT_FIXED} -c ufield3d_lib77.f -o ufield3d_lib77.o

ufield3d_lib.o : ufield3d_lib.f
	${FC} ${OPTS} ${FORMAT_FREE} -c ufield3d_lib.f -o ufield3d_lib.o

ufield3d_class.o : ufield3d_class.f spect3d_class.o perrors_class.o \
                   parallel_pipe_class.o ufield3d_lib.o hdf5io_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c ufield3d_class.f -o ufield3d_class.o

part2d_class.o : part2d_class.f fdist2d_class.o ufield2d_class.o spect2d_class.o \
                  perrors_class.o parallel_pipe_class.o part2d_lib.o hdf5io_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c part2d_class.f -o part2d_class.o

part2d_lib77.o : part2d_lib77.f
	${FC} ${OPTS} ${FORMAT_FIXED} -c part2d_lib77.f -o part2d_lib77.o

part2d_lib.o : part2d_lib.f
	${FC} ${OPTS} ${FORMAT_FREE} -c part2d_lib.f -o part2d_lib.o

fdist2d_class.o : fdist2d_class.f spect2d_class.o perrors_class.o \
                   parallel_pipe_class.o 
	${FC} ${OPTS} ${FORMAT_FREE} -c fdist2d_class.f -o fdist2d_class.o

fpois2d_class.o : fpois2d_class.f ufield2d_class.o spect2d_class.o perrors_class.o \
                   parallel_pipe_class.o fpois2d_lib.o
	${FC} ${OPTS} ${FORMAT_FREE} -c fpois2d_class.f -o fpois2d_class.o

fpois2d_lib77.o : fpois2d_lib77.f
	${FC} ${OPTS} ${FORMAT_FIXED} -c fpois2d_lib77.f -o fpois2d_lib77.o

fpois2d_lib.o : fpois2d_lib.f
	${FC} ${OPTS} ${FORMAT_FREE} -c fpois2d_lib.f -o fpois2d_lib.o


fft2d_class.o : fft2d_class.f ufield2d_class.o spect2d_class.o perrors_class.o \
                   parallel_pipe_class.o fft2d_lib.o
	${FC} ${OPTS} ${FORMAT_FREE} -c fft2d_class.f -o fft2d_class.o

fft2d_lib77.o : fft2d_lib77.f
	${FC} ${OPTS} ${FORMAT_FIXED} -c fft2d_lib77.f -o fft2d_lib77.o

fft2d_lib.o : fft2d_lib.f
	${FC} ${OPTS} ${FORMAT_FREE} -c fft2d_lib.f -o fft2d_lib.o

ufield2d_lib77.o : ufield2d_lib77.f
	${FC} ${OPTS} ${FORMAT_FIXED} -c ufield2d_lib77.f -o ufield2d_lib77.o

ufield2d_lib.o : ufield2d_lib.f
	${FC} ${OPTS} ${FORMAT_FREE} -c ufield2d_lib.f -o ufield2d_lib.o

ufield2d_class.o : ufield2d_class.f spect2d_class.o perrors_class.o \
                   parallel_pipe_class.o ufield2d_lib.o hdf5io_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c ufield2d_class.f -o ufield2d_class.o

hdf5io_class.o : hdf5io_class.f spect3d_class.o perrors_class.o parallel_pipe_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c hdf5io_class.f -o hdf5io_class.o

spect3d_class.o : spect3d_class.f perrors_class.o parallel_pipe_class.o spect2d_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c spect3d_class.f -o spect3d_class.o

spect2d_class.o : spect2d_class.f perrors_class.o parallel_pipe_class.o
	${FC} ${OPTS} ${FORMAT_FREE} -c spect2d_class.f -o spect2d_class.o

dtimer.o : dtimer.c
	gcc -O3 -std=c99 -c dtimer.c

perrors_class.o : perrors_class.f parallel_pipe_class.o dtimer.o
	${FC} ${OPTS} ${FORMAT_FREE} -c perrors_class.f -o perrors_class.o

parallel_pipe_class.o : parallel_pipe_class.f
	${FC} ${OPTS} ${FORMAT_FREE} -c parallel_pipe_class.f -o parallel_pipe_class.o

parallel_class.o : parallel_class.f
	${FC} ${OPTS} ${FORMAT_FREE} -c parallel_class.f -o parallel_class.o

main.o : main.f parallel_pipe_class.o 
	${FC} ${OPTS} ${FORMAT_FREE} -c main.f -o main.o
