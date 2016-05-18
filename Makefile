#
BLAS_ROOT = /opt/OpenBLAS
BLAS_INC_DIR = $(BLAS_ROOT)/include
BLAS_LIB_DIR = $(BLAS_ROOT)/lib
BLAS_LIBS = -lopenblas_seq
#
PLASMA_ROOT = /opt/PLASMA
PLASMA_INC_DIR = $(PLASMA_ROOT)/include
PLASMA_LIB_DIR = $(PLASMA_ROOT)/lib
PLASMA_LIBS = -lplasma -lcoreblas -lquark
#
OMPI_ROOT = /opt/openmpi
OMPI_INC_DIR = $(OMPI_ROOT)/include
OMPI_LIB_DIR = $(OMPI_ROOT)/lib
OMPI_BIN_DIR = $(OMPI_ROOT)/bin
#
TMATRIX_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/TileMatrix
TMATRIX_INC_DIR = $(TMATRIX_ROOT)
TMATRIX_LIB_DIR = $(TMATRIX_ROOT)
TMATRIX_LIBS = -lTileMatrix
#
COREBLAS_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/CoreBlas
COREBLAS_INC_DIR = $(COREBLAS_ROOT)
COREBLAS_LIB_DIR = $(COREBLAS_ROOT)
COREBLAS_LIBS = -lCoreBlasTile
#
PROGRESS_ROOT = /Users/stomo/WorkSpace/TileAlgorithm/ProgressTable
PROGRESS_INC_DIR = $(PROGRESS_ROOT)
PROGRESS_LIB_DIR = $(PROGRESS_ROOT)
PROGRESS_LIBS = -lProgress
#
VT_INC = $(OMPI_INC_DIR)/vampirtrace

# for Debug
#CXX = mpic++
CXX = $(OMPI_BIN_DIR)/mpic++ 
CXXFLAGS =	-DDEBUG -g -fopenmp -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(PROGRESS_INC_DIR) -I$(TMATRIX_INC_DIR) -I$(COREBLAS_INC_DIR)

# for Performance Evaluation
#CXXFLAGS =	-O3 -fopenmp -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(PROGRESS_INC_DIR) -I$(TMATRIX_INC_DIR) -I$(COREBLAS_INC_DIR)

# for Trace
#CXX = mpicxx-vt -vt:inst manual
#CXXFLAGS =	-DVTRACE -fopenmp -I$(BLAS_INC_DIR) -I$(PLASMA_INC_DIR) -I$(PROGRESS_INC_DIR) -I$(TMATRIX_INC_DIR) -I$(COREBLAS_INC_DIR)

OBJS = MPI_Tile.o FB_CAQR.o

TARGET = FB_CAQR

all: $(TARGET)

FB_CAQR : $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) \
				-L$(PROGRESS_LIB_DIR) $(PROGRESS_LIBS) \
				-L$(TMATRIX_LIB_DIR) $(TMATRIX_LIBS) \
				-L$(COREBLAS_LIB_DIR) $(COREBLAS_LIBS) \
				-L$(PLASMA_LIB_DIR) $(PLASMA_LIBS) \
				-L$(BLAS_LIB_DIR) $(BLAS_LIBS)

.cpp.o :
	$(CXX) $(CXXFLAGS) -c $<

clean:
	rm -f Makefile~ *.cpp~ *.h~ *.o $(LLOBJS) $(TARGET)
