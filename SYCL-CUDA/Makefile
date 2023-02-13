CUDA = --expt-extended-lambda -gencode arch=compute_75,code=[sm_75,compute_75] -gencode arch=compute_60,code=[sm_60,compute_60]

CPPFLAGS = -std=c++14 -O2

INCLUDES = -I.

all: builder

builder:
	nvcc ${CPPFLAGS} ${CUDA} ${INCLUDES} -ltbb mainBuilder.cu -o $@.out

builder64:
	nvcc ${CPPFLAGS} ${CUDA} ${INCLUDES} -ltbb mainBuilder.cu -o $@.out

builderDebug:
	nvcc ${CPPFLAGS} ${CUDA} ${INCLUDES} -ltbb -DDEBUG mainBuilder.cu -o $@.out

sycl-cpu: mainBuilder.cpp sycl_octree/octree_traverse.inl basic/traverse_sycl.h
	icpx -fsycl -DCPU -std=c++17 -O2 ${INCLUDES} -o $@ mainBuilder.cpp -tbb -DMEMO -Wno-deprecated-declarations

#without memoization
sycl-cpu-nomemo: mainBuilder.cpp sycl_octree/octree_traverse.inl basic/traverse_sycl.h
	icpx -fsycl -DCPU -std=c++17 -O2 ${INCLUDES} -o $@ mainBuilder.cpp -tbb -Wno-deprecated-declarations

builderSYCLnvidia:
	clang++ -DGPU -DDEVICE -DNVIDIA -std=c++17 -O2 -fsycl -fsycl-targets=nvptx64-nvidia-cuda -Xsycl-target-backend '--expt-extended-lambda' -Xsycl-target-backend '--cuda-gpu-arch=sm_75' -I. mainBuilder.cpp -ltbb -o sycl-builder.out

builderSYCLcuda:
	clang++ -DGPU -DDEVICE -std=c++17 -O2 --cuda-gpu-arch=sm_75 -I/opt/intel/llvm/build/include/sycl -I/users/fmlopez/thrust -I/users/fmlopez/thrust/dependencies/cub ${INCLUDES} mainBuilder.cu -L/opt/intel/llvm/build/lib -lsycl -ltbb -L/usr/local/cuda-11.0/lib64 -lcudart -o sycl-builder.out

example:
	nvcc -ltbb ${CPPFLAGS} ${CUDA} main.cu -o $@.out

example2d:
	nvcc -ltbb ${CPPFLAGS} ${CUDA} main2d.cu -o $@.out

example2dx:
	nvcc -ltbb ${CPPFLAGS} ${CUDA} main2dx.cu -o $@.out

example2ds:
	nvcc -ltbb ${CPPFLAGS} ${CUDA} main2ds.cu -o $@.out

example2dmin:
	nvcc -ltbb ${CPPFLAGS} ${CUDA} main2dmin.cu -o $@.out

example2dprefix:
	nvcc -ltbb ${CPPFLAGS} ${CUDA} main2dprefix.cu -o $@.out

example2dcount:
	nvcc -ltbb ${CPPFLAGS} ${CUDA} main2dcount.cu -o $@.out

examplePrueba:
	nvcc -ltbb ${CPPFLAGS} ${CUDA} mainPrueba.cu -o $@.out

check:
	nvcc -ltbb ${CPPFLAGS} -DCHECK ${CUDA} main.cu -o $@.out

debug:
	nvcc ${CPPFLAGS} -DDEBUG -DCHECK ${CUDA} main.cu -o $@.out

debug2d:
	nvcc -ltbb ${CPPFLAGS} -DDEBUG ${CUDA} main2d.cu -o $@.out

debug2ds:
	nvcc ${CPPFLAGS} -DDEBUG ${CUDA} main2ds.cu -o $@.out

debug2dmin:
	nvcc -ltbb ${CPPFLAGS} -DDEBUG ${CUDA} main2dmin.cu -o $@.out

debug2dprefix:
	nvcc -ltbb ${CPPFLAGS} -DDEBUG ${CUDA} main2dprefix.cu -o $@.out

clean:
	rm -f *.out

# clang++ -std=c++17 --cuda-gpu-arch=sm_75 -I/opt/intel/llvm/build/include/sycl -g vec_add.cu -L/opt/intel/llvm/build/lib -lOpenCL -lsycl -L/usr/local/cuda-11.0/lib64 -lcudart -o vec_add.exe

# clang++ -DUSE_GPU -fsycl -fsycl-targets=nvptx64-nvidia-cuda-sycldevice -Xsycl-target-backend '--cuda-gpu-arch=sm_75' main_prueba.cpp -lsycl -lm -o prueba-out

# Con esto estamos utilizando la Nvidia como un dispositivo OpenCL
# clang++ -DGPU -DDEVICE -fsycl -fsycl-targets=nvptx64-nvidia-cuda-sycldevice -Xsycl-target-backend '--cuda-gpu-arch=sm_75' -I/users/fmlopez/new-lbvh mainBuilder.cpp -lsycl -ltbb -o sycl-builder.out