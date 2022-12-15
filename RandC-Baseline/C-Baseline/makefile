
PRJ_PATH = /home/felipe/Escritorio/Beca_CiTIUS/Proyecto_LiDAR/algoritmo_OWM_LiDAR


INC_DIR_UTILS = utils
SRC_DIR_UTILS = utils


CC = gcc
#CFLAGS =-I/home/remoto/oscar.garcia/datosLidar/LiDAR/gsl/include -L/home/remoto/oscar.garcia/datosLidar/LiDAR/gsl/lib -fopenmp -g -lm -O3 -lgsl -lgslcblas -Wno-unused-result -Wimplicit -lrt
#In Ubuntu -lbsd is needed for strlcpy (#include <bsd/string.h>) sudo apt-get install libbsd-dev
#gsl sudo apt-get install libgsl0-dev
CFLAGS = -fopenmp -g3 -lm -O2 -lgsl -lgslcblas -Wno-unused-result -Wimplicit -lrt -lbsd -llas_c

CFLAGS_CESGA = -fopenmp -g3 -lm -O3 -Wno-unused-result -Wimplicit -lrt -std=c11

INCS_UTILS := $(shell find $(INC_DIR_UTILS) -name '*.h')
SRCS_UTILS := $(shell find $(SRC_DIR_UTILS) -name '*.c')

main : cesga cesga_serie

shared: $(INCS_UTILS) $(SRCS_UTILS)
	$(CC) -c $(SRCS_UTILS) $(INCS_UTILS) $(CFLAGS)

# landing_v8: landing_v8.c landing_v3_plane_gsl.c landing_v3_plane.h landing_approximation.h landing_approximation.c $(INCS_UTILS) $(SRCS_UTILS)
# 	$(CC) landing_v8.c landing_v3_plane_gsl.c landing_v3_plane.h landing_approximation.h landing_approximation.c $(INCS_UTILS) $(SRCS_UTILS) -o landing_v8  $(CFLAGS)

OWM: OWM.c  include/environment.h include/environment.c
	$(CC) OWM.c include/environment.h include/environment.c -o func_OWM.o $(CFLAGS)

datos: datosCesga.c  include/environment.h include/environment.c
	$(CC) datosCesga.c include/environment.h include/environment.c -o datos.o $(CFLAGS)

datos_xyz: datos_xyz.c
	$(CC) datos_xyz.c -o datos_xyz.o

cesga: cesgaOWM.c  include/environment.h include/environment.c
		$(CC) cesgaOWM.c include/environment.h include/environment.c -o parallel.o $(CFLAGS_CESGA)

cesga_serie: cesgaOWM_serie.c  include/environment.h include/environment.c
		$(CC) cesgaOWM_serie.c include/environment.h include/environment.c -o serie.o $(CFLAGS_CESGA)

cesga_out: cesgaOWM_outlayer.c  include/environment.h include/environment.c
		$(CC) cesgaOWM_outlayer.c include/environment.h include/environment.c -o cesga_outlayer.o $(CFLAGS_CESGA)

parallel: OWM_parallel.c  include/environment.h include/environment.c
	$(CC) OWM_parallel.c include/environment.h include/environment.c -o parallel_func_OWM.o $(CFLAGS)

parallelv4: OWM_parallelv4.c include/environment.h include/environment.c
	$(CC) OWM_parallelv4.c include/environment.h include/environment.c -o v4_parallel_func_OWM.o $(CFLAGS)

v6: OWMv6.c  $(INCS_UTILS) $(SRCS_UTILS)
	$(CC) OWMv6.c $(INCS_UTILS) $(SRCS_UTILS) -o func_OWM.o $(CFLAGS)

all: OWM parallel

clean:
	rm -f *.o
