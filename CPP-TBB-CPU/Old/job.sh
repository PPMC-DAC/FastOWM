#!/bin/bash


# NUM_CORES=1
# REPS=3
# for maxpoints in 8 16 32 64; do
#     ./bin/parallelgpu ./data/INAER_2011_Alcoy_Core 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done

# for maxpoints in 16 32 64; do
#     ./bin/parallelgpu ./data/BABCOCK_2017_Arzua_3B 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done

# for maxpoints in 8 16 32 64; do
#     ./bin/parallelgpu ./data/V21_group1_densified_point_cloud 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done

# for maxpoints in 8 16 32 64; do
#     ./bin/parallelgpu ./data/V19_group1_densified_point_cloud 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done



# NUM_CORES=8
# REPS=4
# for maxpoints in 32 40 48 56 64; do
#     ./bin/parallelgpu ./data/INAER_2011_Alcoy_Core 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done

# for maxpoints in 32 40 48 56 64; do
#     ./bin/parallelgpu ./data/BABCOCK_2017_Arzua_3B 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done

# for maxpoints in 32 40 48 56 64; do
#     ./bin/parallelgpu ./data/V21_group1_densified_point_cloud 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done

# for maxpoints in 32 40 48 56 64; do
#     ./bin/parallelgpu ./data/V19_group1_densified_point_cloud 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done



# NUM_CORES=1
# REPS=3
# for maxpoints in 8 16 32 64; do
#     ./bin/parallelgpu ./data/INAER_2011_Alcoy_Core 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done

# for maxpoints in 16 32 64; do
#     ./bin/parallelgpu ./data/BABCOCK_2017_Arzua_3B 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done

# for maxpoints in 8 16 32 64; do
#     ./bin/parallelgpu ./data/V21_group1_densified_point_cloud 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done

# for maxpoints in 8 16 32 64; do
#     ./bin/parallelgpu ./data/V19_group1_densified_point_cloud 10 20 0.8 $NUM_CORES $REPS 0 0 $maxpoints
# done



# NUM_CORES=8
# REPS=1
# for maxpoints in 16 32 48 64; do
#     ./bin/parallelgpu -i ./data/INAER_2011_Alcoy_Core -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r 0 -b 0.2 -s $maxpoints
# done

# for maxpoints in 16 32 48 64; do
#     ./bin/parallelgpu -i ./data/BABCOCK_2017_Arzua_3B -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r 0 -b 0.2 -s $maxpoints
# done

# for maxpoints in 16 32 48 64; do
#     ./bin/parallelgpu -i ./data/V21_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r 0 -b 0.2 -s $maxpoints
# done

# for maxpoints in 16 32 48 64; do
#     ./bin/parallelgpu -i ./data/V19_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r 0 -b 0.2 -s $maxpoints
# done



# NUM_CORES=8
# REPS=2
# for radius in 0.1 0.2 0.3 0.4; do
#     ./bin/parallelgpu -i ./data/INAER_2011_Alcoy_Core -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r $radius -b 0.05 -s 0
# done

# for radius in 0.1 0.2 0.3 0.4; do
#     ./bin/parallelgpu -i ./data/BABCOCK_2017_Arzua_3B -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r $radius -b 0.05 -s 0
# done

# for radius in 0.1 0.2 0.3 0.4; do
#     ./bin/parallelgpu -i ./data/V21_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r $radius -b 0.05 -s 0
# done

# for radius in 0.1 0.2 0.3 0.4; do
#     ./bin/parallelgpu -i ./data/V19_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r $radius -b 0.05 -s 0
# done


# NUM_CORES=8
# REPS=20
# for maxsize in 32 48 64 128; do
#         ./bin/parallelgpu -i ./data/INAER_2011_Alcoy_Core -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r 0 -b 0.05 -s $maxsize -L 2
# done

# for maxsize in 32 48 64 128; do
#         ./bin/parallelgpu -i ./data/BABCOCK_2017_Arzua_3B -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r 0 -b 0.05 -s $maxsize -L 2
# done

# for maxsize in 32 48 64 128; do
#         ./bin/parallelgpu -i ./data/V21_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r 0 -b 0.05 -s $maxsize -L 2
# done

# for maxsize in 32 48 64 128; do
#         ./bin/parallelgpu -i ./data/V19_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n $NUM_CORES -l $REPS -r 0 -b 0.05 -s $maxsize -L 2
# done



# for maxpoints in 256 512 1024 2048; do
#         # for limit in 32 64 128; do
#                 for nlevel in 5 6 7 8; do
#                         ./bin/parallelgpu -i ./data/INAER_2011_Alcoy_Core -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxpoints -L $nlevel -l 3
#                 done
#         # done
# done
# for maxpoints in 256 512 1024 2048; do
#         # for limit in 32 64 128; do
#                 for nlevel in 5 6 7 8; do
#                         ./bin/parallelgpu -i ./data/BABCOCK_2017_Arzua_3B -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxpoints -L $nlevel -l 3
#                 done
#         # done
# done
# for maxpoints in 256 512 1024 2048; do
#         # for limit in 32 64 128; do
#                 for nlevel in 5 6 7 8; do
#                         ./bin/parallelgpu -i ./data/V21_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxpoints -L $nlevel -l 3
#                 done
#         # done
# done
# for maxpoints in 256 512 1024 2048; do
#         # for limit in 32 64 128; do
#                 for nlevel in 5 6 7 8; do
#                         ./bin/parallelgpu -i ./data/V19_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxpoints -L $nlevel -l 3
#                 done
#         # done
# done


# for radius in 0.2 0.5 1.0 2.0; do
#         # for limit in 32 64 128; do
#                 for nlevel in 2 3 4 5 6 7; do
#                         ./bin/parallelgpu -i ./data/INAER_2011_Alcoy_Core -W 10 -B 20 -O 0.8 -n 9 -r $radius -b 0 -s 0 -L $nlevel -l 5
#                 done
#         # done
# done
# for radius in 0.2 0.5 1.0 2.0; do
#         # for limit in 32 64 128; do
#                 for nlevel in 4 5 6 7; do
#                         ./bin/parallelgpu -i ./data/BABCOCK_2017_Arzua_3B -W 10 -B 20 -O 0.8 -n 9 -r $radius -b 0 -s 0 -L $nlevel -l 5
#                 done
#         # done
# done
# for radius in 0.1 0.2 0.5 1.0; do
#         # for limit in 32 64 128; do
#                 for nlevel in 2 3 4 5 6; do
#                         ./bin/parallelgpu -i ./data/V21_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n 9 -r $radius -b 0 -s 0 -L $nlevel -l 5
#                 done
#         # done
# done
# for radius in 0.1 0.2 0.5 1.0; do
#         # for limit in 32 64 128; do
#                 for nlevel in 2 3 4 5 6; do
#                         ./bin/parallelgpu -i ./data/V19_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n 9 -r $radius -b 0 -s 0 -L $nlevel -l 5
#                 done
#         # done
# done


# for maxsize in 8 16 32; do
#         for chunk in 32768 65536 131072 262144 ; do
#                 # for nlevel in 4 5 6 7; do
#                 ./bin/paralleldynidx -i ./data/INAER_2011_Alcoy_Core -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxsize -L 7 -l 3 -c $chunk
#         done
# done

# for maxsize in 8 16 32; do
#         for chunk in 32768 65536 131072 262144 ; do
#                 # for nlevel in 4 5 6 7; do
#                 ./bin/paralleldynidx -i ./data/BABCOCK_2017_Arzua_3B -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxsize -L 7 -l 3 -c $chunk
#         done
# done

# for maxsize in 8 16 32; do
#         for chunk in 32768 65536 131072 262144 ; do
#                 # for nlevel in 4 5 6 7; do
#                 ./bin/paralleldynidx -i ./data/V21_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxsize -L 7 -l 3 -c $chunk
#         done
# done

# for maxsize in 8 16 32; do
#         for chunk in 32768 65536 131072 262144 ; do
#                 # for nlevel in 4 5 6 7; do
#                 ./bin/paralleldynidx -i ./data/V19_group1_densified_point_cloud -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxsize -L 7 -l 3 -c $chunk
#         done
# done


# for maxsize in 64 128 256 512 1024; do
#         ./bin/parallelgpu -i ./data/INAER_2011_Alcoy_CoreX4 -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxsize -L 7 -l 5
# done
# for maxsize in 64 128 256 512 1024; do
#         ./bin/parallelgpu -i ./data/INAER_2011_Alcoy_CoreX6 -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxsize -L 7 -l 5
# done
for maxsize in 8 16 32 64; do
        ./bin/parallelgpu -i ./data/INAER_2011_Alcoy -W 10 -B 20 -O 0.8 -n 8 -r 0 -b 0 -s $maxsize -L 7 -l 5
done

