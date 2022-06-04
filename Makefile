# Location of the CUDA Toolkit on my machine
CUDA_PATH = /usr/local/cuda

# Location of VTK stuff on my machine
VTK_PATH = /usr/local/include/vtk-9.1/
VTK_LIB_PATH = /usr/local/src/VTK/build/lib

# Compilers
NVCC = $(CUDA_PATH)/bin/nvcc
CC = g++

# This creates binaries targeting the K80
NVCCFLAGS = -arch=sm_37 -c

# Flags for ordinary g++
CFLAGS = -c

# Linker flags
LDFLAGS = 
#-Wl,--copy-dt-needed-entries

# Include path
INCLUDE_PATH  := -I$(CUDA_PATH)/include -I$(VTK_PATH) -I. 

# Libs
LCUDA = -L$(CUDA_PATH)/lib64 -lcuda -lcudart
LNUM = -lm
LGL = -lGL -lGLU -lglut
# Normally, one is supposed to make a VTK app using CMake.  I
# didn't want to fiddle around with CMake building a CUDA app, so
# I instead used regular-old GNU make to build.  This necessates
# that I include the following long list of libs into the build.
LVTK = -L$(VTK_LIB_PATH) \
	-lvtkCommonSystem-9.1 \
	-lvtkCommonCore-9.1 \
	-lvtkCommonColor-9.1 \
	-lvtkCommonMisc-9.1 \
	-lvtkCommonDataModel-9.1 \
	-lvtkCommonExecutionModel-9.1 \
	-lvtkChartsCore-9.1 \
	-lvtkRenderingCore-9.1 \
	-lvtkViewsCore-9.1 \
	-lvtkIOCore-9.1 \
	-lvtkInfovisCore-9.1 \
	-lvtkImagingColor-9.1 \
	-lvtkViewsContext2D-9.1 \
	-lvtkRenderingContext2D-9.1 \
	-lvtkRenderingContextOpenGL2-9.1 \
	-lvtkRenderingGL2PSOpenGL2-9.1 \
	-lvtkRenderingOpenGL2-9.1 \
	-lvtkInteractionStyle-9.1 \
	-lvtksys-9.1


LIBRARIES := $(LCUDA) $(LNUM) $(LGL) $(LVTK)

# Stuff to build
MAIN = heateq2D_jacobi
VTKS = make_vtk_plot
SRCS = $(MAIN).cu $(VTKS).c
INCS = $(MAIN).h $(VTKS).h
OBJS = $(MAIN).o $(VTKS).o
EXES = $(MAIN)

#----------------------------------------------------------
all: heateq2D_jacobi

heateq2D_jacobi: heateq2D_jacobi.o make_vtk_plot.o
	echo "--->  Linking heateq2D_jacobi"
	$(NVCC) $(LIBRARIES) $(LDFLAGS) heateq2D_jacobi.o make_vtk_plot.o  -o heateq2D_jacobi

heateq2D_jacobi.o: *.cu *.h
	echo "--->  Building heateq2D_jacobi.o"
	$(NVCC) $(INCLUDE_PATH) $(NVCCFLAGS) heateq2D_jacobi.cu

make_vtk_plot.o: make_vtk_plot.c make_vtk_plot.h
	echo "--->  Building make_vtk_plot.o"
	$(CC) $(INCLUDE_PATH) $(CFLAGS) make_vtk_plot.c

clean:
	rm -f $(EXES) octave-workspace *.o *~
