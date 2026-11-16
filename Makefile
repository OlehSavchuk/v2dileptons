CXX      ?= g++
CXXFLAGS ?= -O2 -std=c++14 -fopenmp -Wall

SRC = dilepton_yield.cpp dileptons.cpp spectral_function.cpp
HDR = dileptons.h spectral_function.h

# the program (run from this folder, it needs ./TriangFiles/)
dilepton_v2: $(SRC) $(HDR)
	$(CXX) $(CXXFLAGS) $(SRC) -o $@

# compares the per-cell fast path with the original emission rate and times both (needs the ImDvec tables)
tests/test_fast: tests/test_fast.cpp dileptons.cpp spectral_function.cpp $(HDR)
	$(CXX) $(CXXFLAGS) tests/test_fast.cpp dileptons.cpp spectral_function.cpp -o $@

test: tests/test_fast
	./tests/test_fast | grep -v ImDvec

# checks TriangFiles/triangles.dat
check-triangulation:
	python3 tools/triangulation.py

clean:
	rm -f dilepton_v2 tests/test_fast

.PHONY: test check-triangulation clean
