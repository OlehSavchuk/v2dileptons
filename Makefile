CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++14 -fopenmp -Wall

SRC = main.cpp dileptonsclass.cpp spectral_func.cpp

a.out: $(SRC) dileptonsclass.h spectral_func.h
	$(CXX) $(CXXFLAGS) $(SRC) -o $@

clean:
	rm -f a.out
.PHONY: clean

tests/test_fast: tests/test_fast.cpp dileptonsclass.cpp spectral_func.cpp dileptonsclass.h spectral_func.h
	$(CXX) $(CXXFLAGS) tests/test_fast.cpp dileptonsclass.cpp spectral_func.cpp -o $@

# needs the TriangFiles/ImDvec tables; compares the fast per-cell path with the original emission rate
test: tests/test_fast
	./tests/test_fast | grep -v ImDvec
.PHONY: test
