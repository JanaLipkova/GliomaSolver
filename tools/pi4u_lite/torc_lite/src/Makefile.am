# include common defs
torclibdir = $(libdir)
torclib_LIBRARIES = libtorc.a

INC = ../include
AM_CFLAGS = @DEBUG_FLAG@ -DMAX_NVPS=@NVPS@ -DMAX_NODES=@NNODES@ -I../include -I../uth -I.
#-DDBG

libtorc_a_SOURCES = torc_runtime.c torc_queue.c torc_thread.c torc_comm.c torc_server.c torc.c

torcincdir = $(includedir)
#torcinc_HEADERS = $(INC)/torc_data.h $(INC)/torc_internal.h $(INC)/torc.h $(INC)/torc_mpi_internal.h $(INC)/torc_queue.h $(INC)/utils.h $(INC)/ps_config.h $(INC)/locks.h $(INC)/queues.h $(INC)/torcf.h
torcinc_HEADERS = $(INC)/torc.h $(INC)/torcf.h 
#$(INC)/torc_config.h
