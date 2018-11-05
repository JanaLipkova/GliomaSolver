SHELL := /bin/bash

CC = g++
OPENMP_FLAG= -fopenmp
		
CPPFLAGS+= -I..

ifeq "$(config)" "production"
CPPFLAGS += -DNDEBUG
endif

##################	


VPATH := ../source/  
.DEFAULT_GOAL := likelihood


OBJECTS = Likelihood_main.o HGG_Likelihood.o

all: likelihood

likelihood: $(OBJECTS)
	$(CC) $^ -o $@ 

%.o: %.cpp
	$(CC) $(CPPFLAGS) -c $^ -o $@

clean:
	rm -f *.o likelihood

