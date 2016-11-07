# Makefile for kmc analysis C++ extension
# Author : shaozhengjiang<shaozhengjiang@gmail.com>
# Date   : 2015-12-29
# Update : 2016-11-07

CXX = g++
SWIG = swig
SWIGFLAGS = -python -c++ -Wall -Werror
CPPFLAGS = -O3 -fPIC -std=c++0x -c -Wall -Werror
INC = -I/usr/include/python2.7/
LDFLAGS = -shared
EXTNAME = kmc_functions
EXTDIR = ../kynetix/solvers/kmc_plugins/plugin_backends

_kmc_functions.so: kmc_functions.o kmc_functions_wrap.o
	$(CXX) $(LDFLAGS) $^ -o $@
	mkdir build
	mv *.o build

kmc_functions.o: src/kmc_functions.cpp
	$(CXX) $(CPPFLAGS) -I include $< -o $@

kmc_functions_wrap.o: wrap/kmc_functions_wrap.cxx
	$(CXX) $(CPPFLAGS) $(INC) -I include $< -o $@

wrap/kmc_functions_wrap.cxx: wrap/kmc_functions.i
	$(SWIG) $(SWIGFLAGS) -I./include $^
	mv wrap/*py ./

install:
	cp _$(EXTNAME).so $(EXTNAME).py $(EXTDIR)

clean:
	rm -rf build wrap/$(EXTNAME)_wrap.cxx $(EXTNAME).py *.so

uninstall:
	rm -f $(EXTDIR)/$(EXTNAME).py $(EXTDIR)/_$(EXTNAME).so

