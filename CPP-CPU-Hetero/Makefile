UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
CC = gcc-12
CFLAGS_BASELINE = -fopenmp -lm -O3 -std=c11
endif

ifeq ($(UNAME), Linux)
CC = gcc
CFLAGS_BASELINE = -fopenmp -lm -O3 -lrt -std=c11
endif

BIN=bin
SRC=src
APP=application

#CFLAGS =-I/home/remoto/oscar.garcia/datosLidar/LiDAR/gsl/include -L/home/remoto/oscar.garcia/datosLidar/LiDAR/gsl/lib -fopenmp -g -lm -O3 -lgsl -lgslcblas -Wno-unused-result -Wimplicit -lrt
#In Ubuntu -lbsd is needed for strlcpy (#include <bsd/string.h>) sudo apt-get install libbsd-dev
#gsl sudo apt-get install libgsl0-dev
CFLAGS = -fopenmp -g3 -lm -O2 -lgsl -lgslcblas -Wno-unused-result -Wimplicit -lrt -lbsd -llas_c

CPPFLAGS_CESGA = -O3 -fopenmp -lm -ltbb -std=c++11

CFLAGS_DBG = -fopenmp -g3 -lm -O2 -Wall -lrt -std=c11 -qopt-report=5 -qopenmp

PROG=baseline
LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)

$(BIN)/%: $(APP)/OWM_baseline.c $(SRC)/OWM_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE)

$(BIN)/par_baseline: $(APP)/OWM_baseline.c $(SRC)/OWM_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL

bin/sequential: application/cesgaOWM_seq.c src/OWM_functions.c
		@mkdir -p "$(@D)"
		$(CC) application/cesgaOWM_seq.c src/OWM_functions.c -o $@ $(CFLAGS_BASELINE)

bin/qold:
		@mkdir -p "$(@D)"
		g++ application/owm_qold.cpp src/envi_qold.cpp -o $@ -O3 -fopenmp -lm -std=c++11

bin/qmin:
		@mkdir -p "$(@D)"
		g++ application/owm_qmin.cpp src/envi_qmin.cpp -o $@ -O3 -fopenmp -lm -std=c++11

bin/parent:
		@mkdir -p "$(@D)"
		g++ application/owm_parent.cpp src/envi_parent.cpp -o $@ -O3 -fopenmp -lm -std=c++11

bin/unique:
		@mkdir -p "$(@D)"
		g++ application/owm_up.cpp src/envi_up.cpp -o $@ -O3 -fopenmp -lm -std=c++11

bin/parallelc: $(APP)/cesgaOWM.c $(SRC)/OWM_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE)

bin/parallelcpp:
		@mkdir -p "$(@D)"
		g++ application/owm_tfm.cpp src/envi.cpp $(CPPFLAGS_CESGA) -o $@

bin/paralleldpcpp:
		@mkdir -p "$(@D)"
		dpcpp application/owm_gpu.cpp src/envi_gpu.cpp -lOpenCL -lsycl -O3 -lm -ltbb -std=c++11 -o $@

bin/parallelgpu:
		@mkdir -p "$(@D)"
		dpcpp application/owm_gpu.cpp -lsycl -O3 -lm -ltbb -ltbbmalloc_proxy -o $@

bin/paralleldyn:
		@mkdir -p "$(@D)"
		dpcpp application/owm_gpu.cpp -lsycl -O3 -lm -ltbb -ltbbmalloc_proxy -std=c++11 -DDYNAMIC -o $@

bin/paralleldynidx:
		@mkdir -p "$(@D)"
		dpcpp application/owm_gpu.cpp -lsycl -O3 -lm -ltbb -ltbbmalloc_proxy -std=c++11 -DDYNAMIC -DINDEX -o $@

bin/parallelindex:
		@mkdir -p "$(@D)"
		dpcpp application/owm_gpu.cpp -lsycl -O3 -lm -ltbb -ltbbmalloc_proxy -std=c++11 -DINDEX -o $@

bin/stats:
		@mkdir -p "$(@D)"
		dpcpp application/owm_stats.cpp -lOpenCL -lsycl -O3 -lm -ltbb -std=c++11 -o $@

bin/cudaParallel:
		@mkdir -p "$(@D)"
		clang++ application/owm_gpu.cpp -fsycl -fsycl-targets=nvptx64-nvidia-cuda-sycldevice -O3 -lm -ltbb -o $@

clean:
	rm -rf bin/
