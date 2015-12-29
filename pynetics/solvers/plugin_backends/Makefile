# Makefile for kmc analysis C extension
# Author: shaozhengjiang<shaozhengjiang@gmail.com>
# Date : 2015.12.29

CC = gcc
SWIG = swig
SWIGFLAGS = -python 
CCFLAGS = -O3 -fPIC -c 
LDFLAGS = -shared
OBJS = kmc_functions.o kmc_functions_wrap.o

_kmc_functions.so: kmc_functions.o kmc_functions_wrap.o
	$(CC) $(LDFLAGS) $^ -o $@

kmc_functions.o: kmc_functions.c
	$(CC) $(CCFLAGS) $^

kmc_functions_wrap.o: kmc_functions_wrap.c
	$(CC) $(CCFLAGS) $^

kmc_functions_wrap.c: kmc_functions.i
	$(SWIG) $(SWIGFLAGS) $^

clean:
	rm -rf kmc_functions_wrap.c kmc_functions.py _kmc_functions.so *.o *.pyc