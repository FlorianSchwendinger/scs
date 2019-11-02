ifeq ($(OS),Windows_NT)
UNAME = CYGWINorMINGWorMSYS
else
UNAME = $(shell uname -s)
endif

#CC = gcc
# For cross-compiling with mingw use these.
#CC = i686-w64-mingw32-gcc -m32
#CC = x86_64-w64-mingw32-gcc-4.8
CUCC = $(CC) #Don't need to use nvcc, since using cuda blas APIs

# For GPU must add cuda libs to path, e.g.
# export DYLD_LIBRARY_PATH=/usr/local/cuda/lib:$DYLD_LIBRARY_PATH

ifneq (, $(findstring CYGWIN, $(UNAME)))
ISWINDOWS := 1
else
ifneq (, $(findstring MINGW, $(UNAME)))
ISWINDOWS := 1
else
ifneq (, $(findstring MSYS, $(UNAME)))
ISWINDOWS := 1
else
ifneq (, $(findstring mingw, $(CC)))
ISWINDOWS := 1
else
ISWINDOWS := 0
endif
endif
endif
endif

ifeq ($(UNAME), Darwin)
# we're on apple, no need to link rt library
LDFLAGS += -lm
SHARED = dylib
SONAME = -install_name
CULDFLAGS = -L/usr/local/cuda/lib
else
ifeq ($(ISWINDOWS), 1)
# we're on windows (cygwin or msys)
LDFLAGS += -lm
SHARED = dll
SONAME = -soname
#TODO: probably doesn't work:
CULDFLAGS = -L/usr/local/cuda/lib64
else
# we're on a linux system, use accurate timer provided by clock_gettime()
LDFLAGS += -lm -lrt
SHARED = so
SONAME = -soname
CULDFLAGS = -L/usr/local/cuda/lib64
endif
endif

# Add on default CFLAGS
OPT = -O3
override CFLAGS += -g -Wall -pedantic -funroll-loops -I. -Iinclude -Ilinsys $(OPT)
ifneq ($(ISWINDOWS), 1)
override CFLAGS += -fPIC
endif

CULDFLAGS += -lcudart -lcublas -lcusparse
CUDAFLAGS = $(CFLAGS) -I/usr/local/cuda/include -Wno-c++11-long-long # turn off annoying long-long warnings in cuda header files

LINSYS = linsys
DIRSRC = $(LINSYS)/direct
DIRSRCEXT = $(DIRSRC)/external
INDIRSRC = $(LINSYS)/indirect
GPU = $(LINSYS)/gpu

OUT = out
AR = ar
ARFLAGS = rv
ARCHIVE = $(AR) $(ARFLAGS)
RANLIB = ranlib
INSTALL = install

ifeq ($(PREFIX),)
  PREFIX = /usr/local
endif

OPT_FLAGS =
########### OPTIONAL FLAGS ##########
# these can all be override from the command line
# e.g. make DLONG=1 will override the setting below
DLONG = 0
ifneq ($(DLONG), 0)
OPT_FLAGS += -DDLONG=$(DLONG) # use longs rather than ints
endif
CTRLC = 1
ifneq ($(CTRLC), 0)
OPT_FLAGS += -DCTRLC=$(CTRLC) # graceful interrupts with ctrl-c
endif
SFLOAT = 0
ifneq ($(SFLOAT), 0)
OPT_FLAGS += -DSFLOAT=$(SFLOAT) # use floats rather than doubles
endif
NOVALIDATE = 0
ifneq ($(NOVALIDATE), 0)
OPT_FLAGS += -DNOVALIDATE=$(NOVALIDATE)$ # remove data validation step
endif
NOTIMER = 0
ifneq ($(NOTIMER), 0)
OPT_FLAGS += -DNOTIMER=$(NOTIMER) # no timing, times reported as nan
endif
COPYAMATRIX = 1
ifneq ($(COPYAMATRIX), 0)
OPT_FLAGS += -DCOPYAMATRIX=$(COPYAMATRIX) # if normalize, copy A
endif
TEST_GPU_MAT_MUL = 0
ifneq ($(TEST_GPU_MAT_MUL), 0)
OPT_FLAGS += -DTEST_GPU_MAT_MUL=$(TEST_GPU_MAT_MUL) # tests GPU matrix multiply for correctness
endif
GPU_TRANSPOSE_MAT = 1
ifneq ($(GPU_TRANSPOSE_MAT), 0)
OPT_FLAGS += -DGPU_TRANSPOSE_MAT=$(GPU_TRANSPOSE_MAT) # tranpose A mat in GPU memory
endif

### VERBOSITY LEVELS: 0,1,2
EXTRA_VERBOSE = 0
ifneq ($(EXTRA_VERBOSE), 0)
OPT_FLAGS += -DEXTRA_VERBOSE=$(EXTRA_VERBOSE) # extra verbosity level
endif

############ OPENMP: ############
# set USE_OPENMP = 1 to allow openmp (multi-threaded matrix multiplies):
# set the number of threads to, for example, 4 by entering the command:
# export OMP_NUM_THREADS=4

USE_OPENMP = 0
ifneq ($(USE_OPENMP), 0)
  override CFLAGS += -fopenmp
  LDFLAGS += -lgomp
endif

############ SDPS: BLAS + LAPACK ############
# set USE_LAPACK = 1 below to enable solving SDPs
# NB: point the libraries to the locations where
# you have blas and lapack installed

USE_LAPACK = 1
ifneq ($(USE_LAPACK), 0)
  # edit these for your setup:
  BLASLDFLAGS = -lblas -llapack #-lgfortran
  LDFLAGS += $(BLASLDFLAGS)
  OPT_FLAGS += -DUSE_LAPACK

  BLAS64 = 0
  ifneq ($(BLAS64), 0)
  OPT_FLAGS += -DBLAS64=$(BLAS64) # if blas/lapack lib uses 64 bit ints
  endif

  NOBLASSUFFIX = 0
  ifneq ($(NOBLASSUFFIX), 0)
  OPT_FLAGS += -DNOBLASSUFFIX=$(NOBLASSUFFIX) # hack to strip blas suffix
  endif

  BLASSUFFIX = "_"
  ifneq ($(BLASSUFFIX), "_")
  OPT_FLAGS += -DBLASSUFFIX=$(BLASSUFFIX) # blas suffix (underscore usually)
  endif
endif

MATLAB_MEX_FILE = 0
ifneq ($(MATLAB_MEX_FILE), 0)
OPT_FLAGS += -DMATLAB_MEX_FILE=$(MATLAB_MEX_FILE) # matlab mex
endif
PYTHON = 0
ifneq ($(PYTHON), 0)
OPT_FLAGS += -DPYTHON=$(PYTHON) # python extension
endif
USING_R = 0
ifneq ($(USING_R), 0)
OPT_FLAGS += -DUSING_R=$(USING_R) # R extension
endif

# debug to see var values, e.g. 'make print-OBJECTS' shows OBJECTS value
print-%: ; @echo $*=$($*)

override CFLAGS += $(OPT_FLAGS)
CUDAFLAGS += $(OPT_FLAGS)
