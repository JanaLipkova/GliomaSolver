SHELL := /bin/bash

# Select makefile for your computer
hn ?= $(shell hostname)

ifneq "$(findstring kraken,$(shell hostname))" ""
CPPFLAGS += -D KRAKEN
include make.kraken
endif

ifneq "$(findstring mpp,$(shell hostname))" ""
CPPFLAGS += -D LRZ_CLUSTER
include make.lrz
endif

ifneq "$(findstring Jana,$(shell hostname))" ""
CPPFLAGS += -D JANA
include make.jana
endif

# Compilations options
ifeq "$(helmholtz)" "hypre"
	CPPFLAGS += -DHYPRE
endif

ifeq "$(config)" "production"
CPPFLAGS += -DNDEBUG
CPPFLAGS += -D_FMMSILENT
endif

#MRAG grid parameters
dim ?= 3
bpd ?= 8
#bpd ?= 16
#bs ?= 16
bs=16

CPPFLAGS += -D_DIM=$(dim) -D_BLOCKSIZE_=$(bs) -D_BPD_=$(bpd)  

ifeq "$(dim)" "2"
	bsz ?= 1
	CPPFLAGS += -D_BLOCKSIZE_Z_=$(bsz)
endif

ifeq "$(bpd)" "16"
        ml ?= 4
        CPPFLAGS += -D _MAXLEVEL_=$(ml)
endif

ifeq "$(bpd)" "8"
        ml ?= 3
        CPPFLAGS += -D _MAXLEVEL_=$(ml)
endif
	

# Libraries
LIBS += \
	-L$(TBB_LIB_DIR) \
	-ltbb \
	-ltbbmalloc \
	-L$(VTK_LIB_DIR) \
	-lvtkViews \
	-lvtkInfovis \
	-lvtkWidgets \
	-lvtkHybrid \
	-lvtkRendering \
	-lvtkGraphics \
	-lvtkverdict \
	-lvtkImaging \
	-lvtkftgl \
	-lvtkfreetype \
	-lvtkIO \
	-lvtkFiltering \
	-lvtkCommon \
	-lm \
	-lvtkDICOMParser \
	-lvtkmetaio \
	-lvtksqlite \
	-lvtkpng \
	-lvtktiff \
	-lvtkjpeg \
	-lvtkexpat \
	-lvtksys \
	-lvtkexoIIc \
	-lvtkNetCDF \
	-lvtklibxml2 \
	-lvtkzlib \
	-lpthread \
	-ldl \
	$(OPENMP_FLAG)
	
	
ifeq "$(helmholtz)" "hypre"
	CC=mpic++
	LD=$(CC)
	CPPFLAGS += -I$(HYPRE_INC_DIR)
	LIBS += -L$(HYPRE_LIB_DIR) -lHYPRE 
endif	



CPPFLAGS += -I../Glioma/ -I../MRAG/ -I../Glioma/Operators/ -I../Glioma/Tests/
VPATH := ../Glioma/ ../MRAG/ ../MRAG/MRAGcore ../Glioma/Tests/
.DEFAULT_GOAL := brain

########################################################################	
OBJECTS = \
Glioma_main.o \
Glioma_ReactionDiffusion.o\
MRAGBoundaryBlockInfo.o \
Glioma_BrainDeformationTimeRelaxation.o \
dat2VP.o \
MRAGProfiler.o \
MRAGWavelets_StaticData.o \
Test.o\

ifeq "$(helmholtz)" "hypre"
OBJECTS += Glioma_BrainDeformation.o 
OBJECTS += HelmholtzTest.o 
OBJECTS += HelmholtzTest_MPI.o
endif

ifeq "$(tests)" "on"
CPPFLAGS += -D_TESTS
endif

brain: $(OBJECTS)
	$(CC) $^ -o $@ $(LIBS)
	@echo done 

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o 
	rm -f brain 


