UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
CC = gcc-12
CXX= g++-12
CFLAGS_BASELINE = -fopenmp -lm -O3 -std=c11
endif

ifeq ($(UNAME), Linux)
CC = gcc
CXX = g++
CFLAGS_BASELINE = -fopenmp -lm -O3 -lrt -std=c11
endif

BIN=bin
SRC=src
APP=application

CFLAGS = -fopenmp -g3 -lm -O2 -lgsl -lgslcblas -Wno-unused-result -Wimplicit -lrt -lbsd -llas_c

CPPFLAGS_TBB = -O3 -fopenmp -lm -ltbb -std=c++11

CFLAGS_DBG = -fopenmp -g3 -lm -O2 -Wall -lrt -std=c11 -qopt-report=5 -qopenmp


$(BIN)/o2and3: application/owm_optim2and3.cpp src/optim2and3_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ $(CPPFLAGS_TBB) -o $@

$(BIN)/o2and3maxnum: application/owm_optim2and3.cpp src/optim2and3_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ $(CPPFLAGS_TBB) -o $@ -DMAXNUMBER

$(BIN)/paralleldpcpp: application/owm_gpu.cpp src/envi_gpu.cpp
		@mkdir -p "$(@D)"
		icpx  $^ -lOpenCL -fsycl -O3 -lm -ltbb -std=c++17 -o $@

$(BIN)/parallel: application/owm_gpu.cpp src/envi_gpu.hpp
		@mkdir -p "$(@D)"
		icpx application/owm_gpu.cpp -fsycl -O3 -lm -ltbb -ltbbmalloc_proxy -std=c++17 -o $@ -Wno-deprecated-declarations

$(BIN)/parallel-nomemo: application/owm_gpu.cpp src/envi_gpu.hpp
		@mkdir -p "$(@D)"
		icpx application/owm_gpu.cpp -fsycl -O3 -lm -ltbb -ltbbmalloc_proxy -std=c++17 -o $@ -Wno-deprecated-declarations -DNOMEMO

$(BIN)/parallelmaxnum: application/owm_gpu.cpp src/envi_gpu.hpp
		@mkdir -p "$(@D)"
		icpx application/owm_gpu.cpp -fsycl -O3 -lm -ltbb -ltbbmalloc_proxy -std=c++17 -o $@ -DMAXNUMBER

$(BIN)/paralleldyn: application/owm_gpu.cpp src/envi_gpu.hpp 
		@mkdir -p "$(@D)"
		icpx application/owm_gpu.cpp -fsycl -O3 -lm -ltbb -ltbbmalloc_proxy -std=c++17 -DDYNAMIC -o $@

$(BIN)/paralleldynidx: application/owm_gpu.cpp src/envi_gpu.hpp
		@mkdir -p "$(@D)"
		icpx application/owm_gpu.cpp -fsycl -O3 -lm -ltbb -ltbbmalloc_proxy -std=c++17 -DDYNAMIC -DINDEX -o $@

$(BIN)/parallelindex: application/owm_gpu.cpp src/envi_gpu.hpp
		@mkdir -p "$(@D)"
		icpx application/owm_gpu.cpp -fsycl -O3 -lm -ltbb -ltbbmalloc_proxy -std=c++17 -DINDEX -o $@

$(BIN)/cudaParallel: application/owm_gpu.cpp src/envi_gpu.hpp
		@mkdir -p "$(@D)"
		clang++ application/owm_gpu.cpp -fsycl -fsycl-targets=nvptx64-nvidia-cuda -ltbbmalloc_proxy -O3 -lm -ltbb -o $@

clean:
	rm -rf bin/
