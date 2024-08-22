CXX = icpx
CXXFLAGS = -fsycl -O3 -xCORE-AVX2
TRASH = main ./Result*
.PHONY: clean

all: build run

build:
	$(CXX) $(CXXFLAGS) Main_.cpp
	
run:
	time ./a.out

clean:
	$(RM) -r $(TRASH)


