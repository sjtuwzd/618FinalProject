CC=gcc
DEBUG=0
INSTRUMENT=1
CFLAGS=-g -O3 -Wall -DDEBUG=$(DEBUG) -DTRACK=$(INSTRUMENT) -std=gnu99 -I /home/zhengdaw/gsl/include -L /home/zhengdaw/gsl/lib -lgsl -L /home/zhengdaw/gsl/lib -lgslcblas
CFILES = main.c file_helper.c
HFILES = file_helper.h
LDFLAGS= -lm -lcurl 

all: port_analyzer


port_analyzer: $(CFILES) $(HFILES)
	$(CC) $(CFLAGS) -o port_analyzer $(CFILES) $(LDFLAGS)

clean:
	rm -f port_analyzer
