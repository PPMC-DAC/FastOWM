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

PROG=baseline par_baseline sequential o0parallel o1qtree o2par o3memoA o3memo o3memoB o4minrad o4maxnum stats statsmaxnum rev1 rev1static rev1guided rev1collap rev1cstatic rev1cguided rev1tasks baseline-rev1 baseline-rev1collap
LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)

$(BIN)/%: $(APP)/OWM_baseline.c $(SRC)/OWM_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE)

$(BIN)/par_baseline: $(APP)/OWM_baseline.c $(SRC)/OWM_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL

$(BIN)/sequential: application/cesgaOWM_seq.c src/OWM_functions.c
		@mkdir -p "$(@D)"
		$(CC) application/cesgaOWM_seq.c src/OWM_functions.c -o $@ $(CFLAGS_BASELINE)

$(BIN)/o0parallel: application/owm_optim0.cpp src/optim0_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ $(CPPFLAGS_TBB) -DPARALLEL

$(BIN)/o1qtree: application/owm_optim1_qtree.cpp src/optim1_qtree_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ -O3 -fopenmp -lm -std=c++11 

# scheduling like baseline
$(BIN)/baseline-rev1: application/owm_baseline-rev1.c src/baseline-rev1_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL
# scheduling with collapse
$(BIN)/baseline-rev1collap: application/owm_baseline-rev1.c src/baseline-rev1_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL -DCOLLAPSE

# scheduling like O1
$(BIN)/rev1: application/owm_rev1.cpp src/rev1_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ -O3 -fopenmp -lm -std=c++11 
$(BIN)/rev1static: application/owm_rev1.cpp src/rev1_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ -O3 -fopenmp -lm -std=c++11 -DSTATIC
$(BIN)/rev1guided: application/owm_rev1.cpp src/rev1_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ -O3 -fopenmp -lm -std=c++11 -DGUIDED
# scheduling with collapse
$(BIN)/rev1collap: application/owm_rev1.cpp src/rev1_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ -O3 -fopenmp -lm -std=c++11 -DCOLLAPSE
$(BIN)/rev1cstatic: application/owm_rev1.cpp src/rev1_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ -O3 -fopenmp -lm -std=c++11 -DSTATIC -DCOLLAPSE
$(BIN)/rev1cguided: application/owm_rev1.cpp src/rev1_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ -O3 -fopenmp -lm -std=c++11 -DGUIDED -DCOLLAPSE
# tasks instead of scheduling
$(BIN)/rev1tasks: application/owm_rev1.cpp src/rev1_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ -O3 -fopenmp -lm -std=c++11 -DTASKS

$(BIN)/o2par: application/owm_optim2.cpp src/optim2_func.cpp include/optim2_func.hpp
		@mkdir -p "$(@D)"
		$(CXX) application/owm_optim2.cpp -o $@ -O3 -lm -std=c++17 -ltbb -ltbbmalloc_proxy -lrt
#-DDEBUG -Wall -pedantic -g -Wextra

$(BIN)/o3memoA: application/owm_optim2.cpp src/optim2_func.cpp include/optim2_func.hpp
		@mkdir -p "$(@D)"
		$(CXX) application/owm_optim2.cpp -o $@ -O3 -lm -std=c++17 -ltbb -ltbbmalloc_proxy -lrt -DMEMOA

$(BIN)/o3memo: application/owm_optim3.cpp src/optim3_func.cpp include/optim3_func.hpp
		@mkdir -p "$(@D)"
		$(CXX) application/owm_optim3.cpp -o $@ -O3 -lm -std=c++17 -ltbb -ltbbmalloc_proxy -lrt

$(BIN)/o3memog: application/owm_optim3.cpp src/optim3_func.cpp include/optim3_func.hpp
		@mkdir -p "$(@D)"
		$(CXX) application/owm_optim3.cpp -o $@ -O3 -lm -std=c++17 -ltbb -lrt -g

$(BIN)/o3memoi: application/owm_optim3.cpp src/optim3_func.cpp include/optim3_func.hpp
		@mkdir -p "$(@D)"
		icpx application/owm_optim3.cpp -o $@ -O3 -lm -std=c++17 -ltbb -lrt -g

$(BIN)/o3memoB: application/owm_optim3.cpp src/optim3_func.cpp include/optim3_func.hpp
		@mkdir -p "$(@D)"
		$(CXX) application/owm_optim3.cpp -o $@ -O3 -lm -std=c++17 -ltbb -ltbbmalloc_proxy -lrt -DNOMEMOA

$(BIN)/o4minrad: application/owm_optim3.cpp src/optim3_func.cpp include/optim3_func.hpp
		@mkdir -p "$(@D)"
		$(CXX) application/owm_optim3.cpp -o $@ -O3 -lm -std=c++17 -ltbb -ltbbmalloc_proxy -lrt

$(BIN)/o4maxnum: application/owm_optim3.cpp src/optim3_func.cpp include/optim3_func.hpp
		@mkdir -p "$(@D)"
		$(CXX) application/owm_optim3.cpp -o $@ -O3 -lm -std=c++17 -ltbb -ltbbmalloc_proxy -lrt -DMAXNUMBER
#-pg -g

$(BIN)/stats: application/owm_stats.cpp include/stats_func.hpp
		@mkdir -p "$(@D)"
		$(CXX) application/owm_stats.cpp -O3 -lm -ltbb -std=c++17 -ltbbmalloc_proxy -o $@

$(BIN)/statsmaxnum: application/owm_stats.cpp include/stats_func.hpp
		@mkdir -p "$(@D)"
		$(CXX) application/owm_stats.cpp -O3 -lm -ltbb -std=c++17 -ltbbmalloc_proxy -o $@ -DMAXNUMBER

clean:
	rm -f bin/*
