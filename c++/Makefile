# Makefile for kmc analysis C++ extension
# Author : shaozhengjiang<shaozhengjiang@gmail.com>
# Date   : 2015-12-29
# Update : 2016-06-25

CXX = g++
SWIG = swig
SWIGFLAGS = -python -c++ -Wall -Werror
CPPFLAGS = -O3 -fPIC -std=c++11 -c -Wall -Werror
INC = # -I /usr/local/lib/python2.7/site-packages/numpy/core/include 
LDFLAGS = -shared
EXTNAME = kmc_functions
EXTDIR = ../kynetix/solvers/kmc_plugins/plugin_backends

_kmc_functions.so: build/kmc_functions.o build/kmc_functions_wrap.o
	$(CXX) $(LDFLAGS) $^ -o $@

build/kmc_functions.o: src/kmc_functions.cpp
	$(CXX) $(CPPFLAGS) $(INC) -I include $< -o $@

build/kmc_functions_wrap.o: wrap/kmc_functions_wrap.cxx
	$(CXX) $(CPPFLAGS) $(INC) -I include $< -o $@

wrap/kmc_functions_wrap.cxx: wrap/kmc_functions.i
	$(SWIG) $(SWIGFLAGS) -I./include $^
	mv wrap/*py ./

install:
	cp _$(EXTNAME).so $(EXTNAME).py $(EXTDIR)

clean:
	rm -rf build/* wrap/$(EXTNAME)_wrap.cxx $(EXTNAME).py *.so

uninstall:
	rm -f $(EXTDIR)/$(EXTNAME).py $(EXTDIR)/_$(EXTNAME).so
