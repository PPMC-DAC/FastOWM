#!/bin/bash
# for chunk in 16 32 64 128 256 512; do
#         for threadsXblock in 32 64 96 128 160; do
#                 for blocks in 4 8 16 32; do
#                         ./example2dmin.out data/INAER_2011_Alcoy.xyz $chunk $threadsXblock $blocks
#                 done
#         done
# done
# for chunk in 16 32 64 128 256 512; do
#         for threadsXblock in 32 64 96 128 160; do
#                 for blocks in 4 8 16 32; do
#                         ./example2dmin.out data/INAER_2011_Alcoy_Core.xyz $chunk $threadsXblock $blocks
#                 done
#         done
# done
# for chunk in 16 32 64 128 256 512; do
#         for threadsXblock in 32 64 96 128 160; do
#                 for blocks in 4 8 16 32; do
#                         ./example2dmin.out data/BABCOCK_2017_Arzua_3B.xyz $chunk $threadsXblock $blocks
#                 done
#         done
# done
# for chunk in 16 32 64 128 256 512; do
#         for threadsXblock in 32 64 96 128 160; do
#                 for blocks in 4 8 16 32; do
#                         ./example2dmin.out data/V21_group1_densified_point_cloud.xyz $chunk $threadsXblock $blocks
#                 done
#         done
# done
# for chunk in 16 32 64 128 256 512; do
#         for threadsXblock in 32 64 96 128 160; do
#                 for blocks in 4 8 16 32; do
#                         ./example2dmin.out data/V19_group1_densified_point_cloud.xyz $chunk $threadsXblock $blocks
#                 done
#         done
# done


# for chunk in 4 8 16 32 64; do
#         ./builder.out data/INAER_2011_Alcoy.xyz $chunk
# done

# echo
# echo

# for chunk in 4 8 16 32 64; do
#         ./builder.out data/INAER_2011_Alcoy_Core.xyz $chunk
# done

# echo
# echo

# for chunk in 4 8 16 32 64; do
#         ./builder.out data/BABCOCK_2017_Arzua_3B.xyz $chunk
# done

# echo
# echo

# for chunk in 4 8 16 32 64; do
#         ./builder.out data/V21_group1_densified_point_cloud.xyz $chunk
# done

# echo
# echo

# for chunk in 4 8 16 32 64; do
#         ./builder.out data/V19_group1_densified_point_cloud.xyz $chunk
# done

# echo
# echo

# for chunk in 8 16 32; do
#         ./builder.out data/INAER_2011_Alcoy_CoreX4.xyz
# done


# for threadsXblock in 32 64 96 128 160 192 224 256; do
#         ./example2dmin.out data/BABCOCK_2017_Arzua_3B.xyz 512 $threadsXblock 16
# done
# for blocks in 4 8 16 32; do
#         ./example2dmin.out data/V19_group1_densified_point_cloud.xyz
# done


./builder.out data/INAER_2011_Alcoy.xyz
./builder.out data/INAER_2011_Alcoy_Core.xyz
./builder.out data/BABCOCK_2017_Arzua_3B.xyz
./builder.out data/V21_group1_densified_point_cloud.xyz
./builder.out data/V19_group1_densified_point_cloud.xyz


# for chunk in 4 8 16 32 64; do
#         for factor in 0.5 0.55 0.6 0.65 0.7 0.75 0.8 0.85 0.9 0.95; do
#                 ./builder.out data/INAER_2011_Alcoy.xyz $chunk $factor
#         done
# done

# echo
# echo

# for chunk in 4 8 16 32 64; do
#         for factor in 0.8 0.85 0.9 0.95; do
#                 ./builder.out data/INAER_2011_Alcoy_Core.xyz $chunk $factor
#         done
# done

# echo
# echo

# for chunk in 4 8 16 32 64; do
#         for factor in 0.8 0.85 0.9 0.95; do
#                 ./builder.out data/BABCOCK_2017_Arzua_3B.xyz $chunk $factor
#         done
# done

# echo
# echo

# for chunk in 4 8 16 32 64; do
#         for factor in 0.8 0.85 0.9 0.95; do
#                 ./builder.out data/V21_group1_densified_point_cloud.xyz $chunk $factor
#         done
# done

# echo
# echo

# for chunk in 4 8 16 32 64; do
#         for factor in 0.8 0.85 0.9 0.95; do
#                 ./builder.out data/V19_group1_densified_point_cloud.xyz $chunk $factor
#         done
# done