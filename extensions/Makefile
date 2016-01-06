# Makefile for kmc analysis C extension
# Author: shaozhengjiang<shaozhengjiang@gmail.com>
# Date : 2015.12.29
# ModifyDate: 2016.1.6

CC = gcc
SWIG = swig
SWIGFLAGS = -python 
CCFLAGS = -O3 -fPIC -c # -fopenmp 
LDFLAGS = -shared
EXTNAME = kmc_functions
EXTDIR = ../pynetics/solvers/plugin_backends

_kmc_functions.so: build/kmc_functions.o build/kmc_functions_wrap.o
	$(CC) $(LDFLAGS) $^ -o $@

build/kmc_functions.o: src/kmc_functions.c
	$(CC) $(CCFLAGS) -I include $< -o $@

build/kmc_functions_wrap.o: wrap/kmc_functions_wrap.c
	$(CC) $(CCFLAGS) -I include $< -o $@

wrap/kmc_functions_wrap.c: wrap/kmc_functions.i
	$(SWIG) $(SWIGFLAGS) $^
	mv wrap/*py ./

install:
	cp _$(EXTNAME).so $(EXTNAME).py $(EXTDIR)

clean:
	rm -rf build/* wrap/$(EXTNAME)_wrap.c $(EXTNAME).py *.so

uninstall:
	rm -f $(EXTDIR)/$(EXTNAME).py $(EXTDIR)/_$(EXTNAME).so
