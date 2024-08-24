CXX = g++
CXXFLAGS = -O3 -std=c++11
TRASH = ./Result*
.PHONY: clean

all: build run

build:
	$(CXX) $(CXXFLAGS) Main_.cpp
	
run:
	time ./a.out

clean:
	$(RM) -r $(TRASH)


