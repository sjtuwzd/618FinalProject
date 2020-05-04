MPI=-DMPI
MPICC = mpicc
CC=gcc
DEBUG=0
INSTRUMENT=1
CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG) -DTRACK=$(INSTRUMENT) -std=gnu99 -I /home/zhengdaw/gsl/include -L /home/zhengdaw/gsl/lib -lgsl -L /home/zhengdaw/gsl/lib -lgslcblas
#CFILES = main.c file_helper.c portfolio_helper.c
#HFILES = file_helper.h portfolio_helper.h

CFILES = portfolio.c driver_lib.c portfolio_lib.c
HFILES = driver_lib.h portfolio_lib.h
LDFLAGS= -lm -lcurl 

all: port_analyzer port_analyzer_parallel


port_analyzer: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -o port_analyzer $(CFILES) $(LDFLAGS)


port_analyzer_parallel: $(CFILES) $(HFILES)
	$(MPICC) $(CFLAGS) $(MPI) -o port_analyzer_parallel $(CFILES) $(LDFLAGS)


clean:
	rm -f port_analyzer*
