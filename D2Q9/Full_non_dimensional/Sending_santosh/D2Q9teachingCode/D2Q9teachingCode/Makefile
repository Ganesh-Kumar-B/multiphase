CXX = g++
CXXFLAGS = -O3 -std=c++11
TRASH = ./Result/*.vtk
.PHONY: clean

all: build run

build:
	$(CXX) $(CXXFLAGS) stub.C
	
run:
	time ./a.out

clean:
	$(RM) -r $(TRASH)


