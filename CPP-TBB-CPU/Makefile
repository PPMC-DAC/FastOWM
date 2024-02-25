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

PROG=baseline par_baseline sequential o0parallel o1qtree o2par o3memoA o3memo o3memoB o4maxnum stats statsmaxnum rev1 rev1static rev1guided rev1collap rev1cstatic rev1cguided rev1tasks
LIST=$(addprefix $(BIN)/, $(PROG))

all: $(LIST)

$(BIN)/%: $(APP)/OWM_baseline.c $(SRC)/OWM_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE)

$(BIN)/par_baseline: $(APP)/OWM_baseline.c $(SRC)/OWM_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL

# dynamic scheduling like baseline
$(BIN)/rev1: $(APP)/owm_rev1.c $(SRC)/rev1_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL
$(BIN)/rev1statis: $(APP)/owm_rev1.c $(SRC)/rev1_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL -DSTATIC
$(BIN)/rev1guided: $(APP)/owm_rev1.c $(SRC)/rev1_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL -DGUIDED
# dynamic scheduling with collapse
$(BIN)/rev1collap: $(APP)/owm_rev1.c $(SRC)/rev1_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL -DCOLLAPSE
$(BIN)/rev1cstatic: $(APP)/owm_rev1.c $(SRC)/rev1_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL -DCOLLAPSE -DSTATIC
$(BIN)/rev1cguided: $(APP)/owm_rev1.c $(SRC)/rev1_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL -DCOLLAPSE -DGUIDED
# tasks instead of scheduling
$(BIN)/rev1tasks: $(APP)/owm_rev1.c $(SRC)/rev1_functions.c
		@mkdir -p "$(@D)"
		$(CC) $^ -o $@ $(CFLAGS_BASELINE) -DPARALLEL -DTASKS

$(BIN)/sequential: application/cesgaOWM_seq.c src/OWM_functions.c
		@mkdir -p "$(@D)"
		$(CC) application/cesgaOWM_seq.c src/OWM_functions.c -o $@ $(CFLAGS_BASELINE)

$(BIN)/o0parallel: application/owm_optim0.cpp src/optim0_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ $(CPPFLAGS_TBB) -DPARALLEL

$(BIN)/o1qtree: application/owm_optim1_qtree.cpp src/optim1_qtree_func.cpp
		@mkdir -p "$(@D)"
		$(CXX) $^ -o $@ -O3 -fopenmp -lm -std=c++11 
#-DDEBUG

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

$(BIN)/o3memoB: application/owm_optim3.cpp src/optim3_func.cpp include/optim3_func.hpp
		@mkdir -p "$(@D)"
		$(CXX) application/owm_optim3.cpp -o $@ -O3 -lm -std=c++17 -ltbb -ltbbmalloc_proxy -lrt -DNOMEMOA

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

# run all the rev1 tests
runrev1:
	@echo "Running rev1 guide"
	@bin/rev1guide bin/data/AlcoyH 10 20 0.8 8 2 2
	@echo "Running rev1 static"
	@bin/rev1static bin/data/AlcoyH 10 20 0.8 8 2 2
	@echo "Running rev1 tasks"
	@bin/rev1tasks bin/data/AlcoyH 10 20 0.8 8 2 2
	@echo "Running rev1 collapse"
	@bin/rev1collap bin/data/AlcoyH 10 20 0.8 8 2 2
	@echo "Running rev1"
	@bin/rev1 bin/data/AlcoyH 10 20 0.8 8 2 2