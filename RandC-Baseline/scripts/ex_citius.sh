#!/bin/sh
#PBS -l nodes=1:ppn=64,walltime=0:02:00
#PBS -o citius_salida.out
#PBS -e citius_salida.err
#PBS -m ae -M felipepower7@gmail.com
cd /home/local/felipe.munoz/algoritmo_OWM_LiDAR
./cesga_func_OWM.o ../datos/INAER_2011_Alcoy 12 20 0.8 16 5
