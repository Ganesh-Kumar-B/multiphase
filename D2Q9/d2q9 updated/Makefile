CXX = g++
CXXFLAGS = -O3
TRASH = ./Result*
.PHONY: clean

all: build run

build:
	$(CXX) $(CXXFLAGS) stub.C
	
run:
	time ./a.out

clean:
	$(RM) -r $(TRASH)


