# get CUDA COMPUTE_CAPABILITY for the current GPU
COMPUTE_CAPABILITY:=$(shell /usr/local/cuda/extras/demo_suite/deviceQuery | grep 'CUDA Capability' | awk '{print $$6}' | sed 's/\.//g' | head -n 1)

CUDA = --expt-extended-lambda -gencode arch=compute_${COMPUTE_CAPABILITY},code=[sm_${COMPUTE_CAPABILITY},compute_${COMPUTE_CAPABILITY}] -Xcompiler -fPIE
#-gencode arch=compute_60,code=[sm_60,compute_60]

CPPFLAGS = -std=c++17 -O3 

INCLUDES = -I. -I/usr/include/x86_64-linux-gnu/c++/12

all: owm-sycl-cpu owm-sycl-cpu-nomemo owm-sycl-igpu owm-sycl-igpu-nomemo owm-sycl-dgpu owm-sycl-dgpu-nomemo owm-cuda owm-cuda-grid owm-cuda-nomemo owm-sycl-dgpu-heter


owm-sycl-cpu: mainBuilder.cpp sycl_octree/octree_traverse.inl basic/traverse_sycl.h
	icpx -fsycl -DCPU ${CPPFLAGS} ${INCLUDES} -o $@ mainBuilder.cpp -tbb -Wno-deprecated-declarations
	mv $@ bin

owm-sycl-cpug: mainBuilder.cpp sycl_octree/octree_traverse.inl basic/traverse_sycl.h
	icpx -fsycl -DCPU ${CPPFLAGS} ${INCLUDES} -o $@ mainBuilder.cpp -tbb -Wno-deprecated-declarations -g
	mv $@ bin

#without memoization
owm-sycl-cpu-nomemo: mainBuilder.cpp sycl_octree/octree_traverse.inl basic/traverse_sycl.h
	icpx -fsycl -DCPU -DNOMEMO ${CPPFLAGS} ${INCLUDES} -o $@ mainBuilder.cpp -tbb -Wno-deprecated-declarations
	mv $@ bin

owm-sycl-igpu: mainBuilder.cpp
	icpx -DGPU -DDEVICE ${CPPFLAGS} ${INCLUDES} -fsycl mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
	mv $@ bin

owm-sycl-igpug: mainBuilder.cpp
	icpx -DGPU -DDEVICE ${CPPFLAGS} ${INCLUDES} -fsycl mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations -g
	mv $@ bin

owm-sycl-igpu-debug: mainBuilder.cpp
	icpx -DGPU -DDEVICE -DDEBUG -std=c++17 -fdiagnostics-color=always -g ${INCLUDES} -fsycl mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
	mv $@ bin

owm-sycl-igpu-nomemo: mainBuilder.cpp
	icpx -DGPU -DDEVICE -DNOMEMO ${CPPFLAGS} ${INCLUDES} -fsycl mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
	mv $@ bin

# owm-sycl-dgpu: mainBuilder.cpp
# 	clang++ -DGPU -DDEVICE -DNVIDIA ${CPPFLAGS} ${INCLUDES} -fsycl -fsycl-targets=nvptx64-nvidia-cuda mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
# 	mv $@ bin
# #-Xsycl-target-backend '--expt-extended-lambda' -Xsycl-target-backend '--cuda-gpu-arch=sm_75'
owm-sycl-dgpu: mainBuilder.cpp
	icpx -DGPU -DDEVICE -DNVIDIA ${CPPFLAGS} ${INCLUDES} -fsycl -fsycl-targets=nvptx64-nvidia-cuda mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
	mv $@ bin

owm-sycl-dgpu-heter: mainBuilder.cpp
	icpx -DGPU -DDEVICE -DNVIDIA -DHETER ${CPPFLAGS} ${INCLUDES} -fsycl -fsycl-targets=nvptx64-nvidia-cuda mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
	mv $@ bin

owm-sycl-dgpu-ncu: mainBuilder.cpp
	icpx -DGPU -DDEVICE -DNVIDIA ${CPPFLAGS} ${INCLUDES} -fsycl -fsycl-targets=nvptx64-nvidia-cuda mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations -g
	mv $@ bin

# owm-sycl-dgpu-nomemo: mainBuilder.cpp
# 	clang++ -DGPU -DDEVICE -DNVIDIA -DNOMEMO ${CPPFLAGS} ${INCLUDES} -fsycl -fsycl-targets=nvptx64-nvidia-cuda mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
# 	mv $@ bin
owm-sycl-dgpu-nomemo: mainBuilder.cpp
	icpx -DGPU -DDEVICE -DNVIDIA -DNOMEMO ${CPPFLAGS} ${INCLUDES} -fsycl -fsycl-targets=nvptx64-nvidia-cuda mainBuilder.cpp -ltbb -o $@ -Wno-deprecated-declarations
	mv $@ bin

# owm-sycl-dgpu-interop: mainBuilder.cu
# 	clang++ -DGPU -DDEVICE ${CPPFLAGS} ${CXXFLAGS} --cuda-gpu-arch=sm_70 ${INCLUDES} mainBuilder.cu -fsycl -fsycl-targets=nvptx64-nvidia-cuda -L/usr/local/cuda/lib64 -lcudart -ltbb  -o $@ 
# 	mv $@ bin
# #-L/usr/local/cuda-11.0/lib64 -fsycl-targets=nvptx64-nvidia-cuda -lcudart

CXXFLAGS = -resource-dir=/opt/intel/oneapi/dpcpp-ct/2022.2.1/lib/clang/15.0.0 -Xclang -fcuda-allow-variadic-functions -L/opt/intel/oneapi/compiler/2022.2.1/linux/lib -lOpenCL -fno-sycl-libspirv
#-fsycl-libspirv-path=/opt/intel/oneapi/compiler/2023.1.0/linux/lib/clc/remangled-l64-signed_char.libspirv-nvptx64-nvidia-cuda.bc

owm-sycl-dgpu-interop: mainBuilder.cu
	/opt/intel/oneapi/compiler/2023.1.0/linux/bin-llvm/clang++ -DGPU -DDEVICE ${CPPFLAGS} ${CXXFLAGS} --cuda-gpu-arch=sm_${COMPUTE_CAPABILITY} ${INCLUDES} mainBuilder.cu -fsycl -fsycl-targets=nvptx64-nvidia-cuda -L/usr/local/cuda/lib64 -lcudart_static -allow-unsupported-compiler -ltbb -o $@ 
	mv $@ bin
#-L/usr/local/cuda-11.0/lib64 -fsycl-targets=nvptx64-nvidia-cuda -lcudart

paralleldynidx: owm_gpu.cpp envi_gpu.hpp
	@mkdir -p "$(@D)"
	icpx owm_gpu.cpp -fsycl -O3 -lm -ltbb -ltbbmalloc_proxy -std=c++17 -DDYNAMIC -DINDEX -o $@

BACKEND= -Xsycl-target-backend=nvidia_gpu_sm_89 '-allow-unsupported-compiler'

# XCOMPILER=/opt/intel/oneapi/2024.0/bin/compiler/clang++
# XCOMPILER=/opt/intel/oneapi/compiler/2024.0/bin/compiler/clang++
# XCOMPILER=/opt/intel/oneapi/compiler/2023.0.0/linux/bin-llvm/clang++
# XCOMPILER=/opt/intel/oneapi/compiler/2022.2.1/linux/bin-llvm/clang++
# XCOMPILER=/opt/intel/oneapi/compiler/2023.1.0/linux/bin-llvm/clang++
# XCOMPILER=/opt/intel/oneapi/compiler/2023.2.0/linux/bin-llvm/clang++

owm-hello: hello.cu
	icpx hello.cu -o $@ -std=c++17 --cuda-gpu-arch=sm_89 -fsycl -fsycl-targets=nvptx64-nvidia-cuda -resource-dir=/opt/intel/oneapi/dpcpp-ct/latest/opt/dpct/lib/clang/18 -fsycl-libspirv-path=/opt/intel/oneapi/2024.0/lib/clc/remangled-l64-signed_char.libspirv-nvptx64-nvidia-cuda.bc -L/usr/local/cuda-11/lib64 -lcudart -ldl -lrt

owm-hello-nvcc: hello.cu
	nvcc --compiler-bindir ${XCOMPILER} -allow-unsupported-compiler ${CPPFLAGS} ${CUDA} ${INCLUDES} hello.cu -o $@

owm-cuda: mainBuilder.cu
	nvcc ${CPPFLAGS} ${CUDA} ${INCLUDES} -ltbb mainBuilder.cu -o $@
	mv $@ bin

owm-cudag: mainBuilder.cu
	nvcc ${CPPFLAGS} ${CUDA} ${INCLUDES} -ltbb mainBuilder.cu -o $@ -g
	mv $@ bin

owm-cuda-grid: mainBuilder.cu
	nvcc ${CPPFLAGS} ${CUDA} ${INCLUDES} -ltbb mainBuilder.cu -o $@ -DAUTOGRID
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
	rm -f bin/*

# clang++ -std=c++17 --cuda-gpu-arch=sm_75 -I/opt/intel/llvm/build/include/sycl -g vec_add.cu -L/opt/intel/llvm/build/lib -lOpenCL -lsycl -L/usr/local/cuda-11.0/lib64 -lcudart -o vec_add.exe

# clang++ -DUSE_GPU -fsycl -fsycl-targets=nvptx64-nvidia-cuda-sycldevice -Xsycl-target-backend '--cuda-gpu-arch=sm_75' main_prueba.cpp -lsycl -lm -o prueba-out

# Con esto estamos utilizando la Nvidia como un dispositivo OpenCL
# clang++ -DGPU -DDEVICE -fsycl -fsycl-targets=nvptx64-nvidia-cuda-sycldevice -Xsycl-target-backend '--cuda-gpu-arch=sm_75' -I/users/fmlopez/new-lbvh mainBuilder.cpp -lsycl -ltbb -o sycl-builder.out