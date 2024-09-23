CXX = icpx
CXXFLAGS = -O3 -xCORE-AVX2
TRASH = ./Result*
.PHONY: clean

all: build run

build:
	$(CXX) $(CXXFLAGS) stub.C
	
run:
	time ./a.out

clean:
	$(RM) -r $(TRASH)


