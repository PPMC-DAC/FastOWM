CUDA = --expt-extended-lambda -gencode arch=compute_75,code=[sm_75,compute_75] 
#-gencode arch=compute_60,code=[sm_60,compute_60]

CPPFLAGS = -std=c++17 -O2

INCLUDES = -I.

all: owm-sycl-cpu owm-sycl-cpu-nomemo owm-sycl-igpu owm-sycl-igpu-nomemo owm-sycl-dgpu owm-sycl-dgpu-nomemo owm-cuda


owm-sycl-cpu: mainBuilder.cpp sycl_octree/octree_traverse.inl basic/traverse_sycl.h
	icpx -fsycl -DCPU ${CPPFLAGS} ${INCLUDES} -o $@ mainBuilder.cpp -tbb -Wno-deprecated-declarations
	mv $@ bin

#without memoization
owm-sycl-cpu-nomemo: mainBuilder.cpp sycl_octree/octree_traverse.inl basic/traverse_sycl.h
	icpx -fsycl -DCPU -DNOMEMO ${CPPFLAGS} ${INCLUDES} -o $@ mainBuilder.cpp -tbb -Wno-deprecated-declarations
	mv $@ bin

owm-sycl-igpu: mainBuilder.cpp
	icpx -DGPU -DDEVICE ${CPPFLAGS} ${INCLUDES} -fsycl mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
	mv $@ bin

owm-sycl-igpu-nomemo: mainBuilder.cpp
	icpx -DGPU -DDEVICE -DNOMEMO ${CPPFLAGS} ${INCLUDES} -fsycl mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
	mv $@ bin

owm-sycl-dgpu: mainBuilder.cpp
	clang++ -DGPU -DDEVICE -DNVIDIA ${CPPFLAGS} ${INCLUDES} -fsycl -fsycl-targets=nvptx64-nvidia-cuda mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
	mv $@ bin
#-Xsycl-target-backend '--expt-extended-lambda' -Xsycl-target-backend '--cuda-gpu-arch=sm_75'
owm-sycl-dgpu-nomemo: mainBuilder.cpp
	clang++ -DGPU -DDEVICE -DNVIDIA -DNOMEMO ${CPPFLAGS} ${INCLUDES} -fsycl -fsycl-targets=nvptx64-nvidia-cuda mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
	mv $@ bin

owm-sycl-dgpu-interop: mainBuilder.cu
	clang++ -DGPU -DDEVICE ${CPPFLAGS} --cuda-gpu-arch=sm_75 ${INCLUDES} mainBuilder.cu -fsycl -ltbb -L/usr/local/cuda-11.0/lib64 -lcudart -o $@
	mv $@ bin

owm-cuda: mainBuilder.cu
	nvcc ${CPPFLAGS} ${CUDA} ${INCLUDES} -ltbb mainBuilder.cu -o $@
	mv $@ bin

owm-cuda-nomemo: mainBuilder.cu
	nvcc ${CPPFLAGS} ${CUDA} ${INCLUDES} -ltbb mainBuilder.cu -o $@ -DNOMEMO
	mv $@ bin

owm-cuda-Debug: mainBuilder.cu
	nvcc ${CPPFLAGS} ${CUDA} ${INCLUDES} -ltbb -DDEBUG mainBuilder.cu -o $@.out
	mv $@ bin

#-I/users/fmlopez/thrust -I/users/fmlopez/thrust/dependencies/cub -I/opt/intel/llvm/build/include/sycl -L/opt/intel/llvm/build/lib

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