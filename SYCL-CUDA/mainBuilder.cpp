// #define CL_TARGET_OPENCL_VERSION 300
// #define CUDA_NO_HALF 1

// #include <CL/sycl/backend/cuda.hpp>

// #include <cublas_v1.h>
// #include <cublas_v2.h>
// #include <cuda.h>

#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>

// #include <sycl_ordered/sycl_test.h>
#include <sycl_octree/octree_test.h>

int main(int argc, char* argv[])
{

	std::string inputTXT = (argc > 1)? argv[1] : "data/INAER_2011_Alcoy.xyz";
    uint32_t leafSize = (argc > 2)? uint32_t(atoi(argv[2])) : 16;
	float factor = (argc > 3)? atof(argv[3]) : 0.75;

	// bintree_test(inputTXT, leafSize);
	// bintree_traverse(inputTXT, leafSize);
	// bintree64_test(inputTXT, leafSize);
	// bintree64_traverse(inputTXT, leafSize);
	// ordered_test(inputTXT, leafSize);
	// ordered_traverse(inputTXT, leafSize);
	// octree_test(inputTXT, leafSize);
	octree_traverse(inputTXT, leafSize);
	// octree_traverse_heter(inputTXT, leafSize, factor);
	// sycl_test(inputTXT, leafSize);
	// sycl_traverse(inputTXT, leafSize);
	
	return 0;
}
