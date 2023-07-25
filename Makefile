#   Makefile -- creating a shared C++ library which is imported in MATLAB \
    Copyright (c) 2023 Andreas Michael <andreas.michael@it.uu.se> \
    All rights reserved. Use of this source code is governed by a \
    BSD-style license that can be found in the LICENSE file.
	
INCLUDES= -I$${CONDA_PREFIX}/include -I$${CONDA_PREFIX}/include/eigen3 $$(python3 -m pybind11 --includes)
CFLAGS = -O3 -Wall -Werror -fopenmp -fPIC 
LDFLAGS = -std=c++11 -lspqr -lcholmod -lm -fopenmp -fPIC 
SFILE = bindings
DFILE = src/cppFuncs/
SUFFIX = $$(python3.10-config --extension-suffix)

cppFuncs: solvers.o $(SFILE).o
	g++ -shared solvers.o $(SFILE).o -o cppFuncs$(SUFFIX) $(LDFLAGS) 
solvers.o: $(DFILE)solvers.cpp $(DFILE)funcs.hpp
	g++ $(CFLAGS) $(INCLUDES) -c $(DFILE)solvers.cpp -o solvers.o
bindings.o: $(DFILE)$(SFILE).cpp $(DFILE)funcs.hpp
	g++ $(CFLAGS) $(INCLUDES) -c $(DFILE)$(SFILE).cpp -o $(SFILE).o
clean:
	rm -f cppFuncs* *.o
