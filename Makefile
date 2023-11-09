CXX = g++
CXXFLAGS = -Wall -std=c++11
EIGEN_INCLUDE = eigen-3.4.0 # Replace with the actual path to your Eigen installation

.PHONY: all clean

all: hec

hec: hec.cpp
	$(CXX) $(CXXFLAGS) -I$(EIGEN_INCLUDE) -o $@ $<

clean:
	rm -f hec
