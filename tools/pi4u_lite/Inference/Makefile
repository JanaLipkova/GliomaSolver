hn ?= $(shell hostname)

ifneq "$(findstring mpp, $(shell hostname))" ""
export LANG=C
export LC_ALL=C
GSL_INC_DIR=/lrz/sys/libraries/gsl/2.3/include
GSL_LIB_DIR=/lrz/sys/libraries/gsl/2.3/lib
else
GSL_INC_DIR=${HOME}/usr/gsl/include
GSL_LIB_DIR=${HOME}/usr/gsl/lib
endif

#GSLROOT=${HOME}/usr/gsl
#GSLROOT=/opt/local

CC = mpicc
CPP = mpic++ 
MPIF90 = mpif90

CFLAGS = -O0 -g -ggdb -I../pndl/src
#CFLAGS = -O3 -std=c99
CFLAGS += -Wall #-Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-result -Wall # -w #-Wall

CFLAGS += `torc_cflags`

# GSL
CFLAGS  += -I${GSL_INC_DIR}
LDFLAGS += -L${GSL_LIB_DIR} -lgsl -lgslcblas
#CFLAGS  += -I${GSLROOT}/include
#LDFLAGS += -L${GSLROOT}/lib -lgsl -lgslcblas

# PRIMARY
LDFLAGS += `torc_libs` -lpthread -lm 

TARGETS = display_gen engine_cmaes engine_tmcmc
#TARGETS = display_gen engine_tmcmc engine_tmcmc_sn engine_cmaes engine_tmcmc_glioma
#TARGETS += sa_deriv fd_deriv fd_grad
#TARGETS += engine_ss1 engine_ss2 
TARGETS += propagation_tool

all: $(TARGETS)

engine_tmcmc: engine_tmcmc.c engine_tmcmc.h tmcmc_stats.c tmcmc_db.c gnuplot_i.c auxil.c fitfun_glioma.c
	$(CC) $(CFLAGS) -o engine_tmcmc engine_tmcmc.c tmcmc_stats.c tmcmc_db.c gnuplot_i.c auxil.c $(LDFLAGS)

display_gen: display_gen.c gnuplot_i.c
	$(CC) $(CFLAGS) -o display_gen display_gen.c gnuplot_i.c $(LDFLAGS)

engine_cmaes: engine_cmaes.c cmaes.c fitfun_glioma.c
	$(CC) $(CFLAGS) -o engine_cmaes engine_cmaes.c cmaes.c $(LDFLAGS)

engine_ss1: engine_ss1.cpp ss_db.cpp auxil.o fitfun.c
	$(CPP) $(CFLAGS) -o engine_ss1 engine_ss1.cpp ss_db.cpp auxil.o $(LDFLAGS)

engine_ss2: engine_ss2.cpp ss_db.cpp auxil.o fitfun.c
	$(CPP) $(CFLAGS) -o engine_ss2 engine_ss2.cpp ss_db.cpp auxil.o $(LDFLAGS)

auxil.o: auxil.c
	$(CC) $(CFLAGS) -c auxil.c 

engine_tmcmc_sn: engine_tmcmc_sn.c engine_tmcmc.h tmcmc_stats.c tmcmc_db.c gnuplot_i.c auxil.c fitfun.c
	$(MPIF90) $(CFLAGS) -D_TMCMC_SN_ -o engine_tmcmc_sn engine_tmcmc_sn.c tmcmc_stats.c tmcmc_db.c gnuplot_i.c auxil.c ../pndl/src/libpndl.a $(LDFLAGS)

sa_deriv: sa_deriv.c auxil.o fitfun.c
	$(MPIF90) $(CFLAGS) -o sa_deriv sa_deriv.c auxil.o $(LDFLAGS)

fd_deriv: fd_deriv.c auxil.o fitfun.c
	$(MPIF90) $(CFLAGS) -o fd_deriv fd_deriv.c auxil.o ../pndl/src/libpndl.a $(LDFLAGS)

fd_grad: fd_grad.c auxil.o fitfun.c
	$(MPIF90) $(CFLAGS) -o fd_grad fd_grad.c auxil.o ../pndl/src/libpndl.a $(LDFLAGS)

propagation_tool: propagation_tool.c auxil.c fitfun_glioma.c
	$(CC) $(CFLAGS) -o propagation_tool propagation_tool.c auxil.c $(LDFLAGS)

clean:
	rm -f engine_tmcmc engine_tmcmc_glioma fitfun_glioma display_gen engine_cmaes engine_ss1 engine_ss2 engine_tmcmc_sn sa_deriv fd_deriv fd_grad propagation_tool *.o

clear:
	rm -rf lsf.o* *.dat curgen_db*.txt curres_db*.txt P*S*T* runinfo.txt samples*.txt seeds*.txt Pc.txt actparcmaes.par

backup:
	scp -r $(PWD) chatzidp@bulldozer.ethz.ch:
