#!/bin/bash

# ./sycl-builder-cpu-host.out
# ./sycl-builder-cpu-shared.out
# ./sycl-builder-gpu-host.out
# ./sycl-builder-gpu-shared.out
# ./sycl-builder-gpu-device.out





# for chunk in 0.25 0.3 0.35 0.4 0.45 0.5; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/INAER_2011_Alcoy.xyz 16 $chunk
# done

# echo
# echo

# for chunk in 0.25 0.3 0.35 0.4 0.45 0.5; do
#         ./sycl-builder.out data/INAER_2011_Alcoy_Core.xyz 16 $chunk
# done

# echo
# echo

# for chunk in 0.25 0.3 0.35 0.4 0.45 0.5; do
#         ./sycl-builder.out data/BABCOCK_2017_Arzua_3B.xyz 16 $chunk
# done

# echo
# echo

# for chunk in 0.25 0.3 0.35 0.4 0.45 0.5; do
#         ./sycl-builder.out data/V21_group1_densified_point_cloud.xyz 16 $chunk
# done

# echo
# echo

# for chunk in 0.25 0.3 0.35 0.4 0.45 0.5; do
#         ./sycl-builder.out data/V19_group1_densified_point_cloud.xyz 16 $chunk
# done





# for chunk in 0.6 0.65 0.7 0.75 0.8 0.85 0.9; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/INAER_2011_Alcoy.xyz 16 $chunk
# done

# echo
# echo

# for chunk in 0.6 0.65 0.7 0.75 0.8 0.85 0.9; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/INAER_2011_Alcoy_Core.xyz 16 $chunk
# done

# echo
# echo

# for chunk in 0.6 0.65 0.7 0.75 0.8 0.85 0.9; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/BABCOCK_2017_Arzua_3B.xyz 16 $chunk
# done

# echo
# echo

# for chunk in 0.6 0.65 0.7 0.75 0.8 0.85 0.9; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/V21_group1_densified_point_cloud.xyz 16 $chunk
# done

# echo
# echo

# for chunk in 0.6 0.65 0.7 0.75 0.8 0.85 0.9; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/V19_group1_densified_point_cloud.xyz 16 $chunk
# done


# for chunk in 8 16 32 64 128 256; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/INAER_2011_Alcoy.xyz $chunk
# done

# echo
# echo

# for chunk in 8 16 32 64 128 256; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/INAER_2011_Alcoy_Core.xyz $chunk
# done

# echo
# echo

# for chunk in 8 16 32 64 128 256; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/BABCOCK_2017_Arzua_3B.xyz $chunk
# done

# echo
# echo

# for chunk in 8 16 32 64 128 256; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/V21_group1_densified_point_cloud.xyz $chunk
# done

# echo
# echo

# for chunk in 8 16 32 64 128 256; do
#         SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/V19_group1_densified_point_cloud.xyz $chunk
# done




# ./sycl-builder.out data/INAER_2011_Alcoy.xyz
# ./sycl-builder.out data/INAER_2011_Alcoy_Core.xyz
# ./sycl-builder.out data/BABCOCK_2017_Arzua_3B.xyz
# ./sycl-builder.out data/V21_group1_densified_point_cloud.xyz
# ./sycl-builder.out data/V19_group1_densified_point_cloud.xyz

SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/INAER_2011_Alcoy.xyz
SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/INAER_2011_Alcoy_Core.xyz
SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/BABCOCK_2017_Arzua_3B.xyz
SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/V21_group1_densified_point_cloud.xyz
SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/V19_group1_densified_point_cloud.xyz

# SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/INAER_2011_Alcoy_CoreX4.xyz
# SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/INAER_2011_Alcoy_CoreX6.xyz
# SYCL_DEVICE_FILTER=PI_CUDA ./sycl-builder.out data/INAER_2011_Alcoy_CoreX8.xyz
