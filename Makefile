CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++14 -fopenmp -Wall

SRC = main.cpp dileptonsclass.cpp spectral_func.cpp

a.out: $(SRC) dileptonsclass.h spectral_func.h
	$(CXX) $(CXXFLAGS) $(SRC) -o $@

clean:
	rm -f a.out
.PHONY: clean
