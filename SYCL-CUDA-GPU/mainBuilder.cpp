
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>

// #include <sycl_ordered/sycl_test.h>
#include <sycl_octree/octree_test.h>

int main(int argc, char* argv[])
{
	std::string inputTXT = (argc > 1)? argv[1] : "data/INAER_2011_Alcoy.xyz";
    uint32_t maxNumber = (argc > 2)? uint32_t(atoi(argv[2])) : 16; //Number of LiDAR points per leaf-node
	float factor = (argc > 3)? atof(argv[3]) : 0.75; //Only used for octree_traverse_heter that traverses the tree in CPU+GPU
	// OWM window size
	uint32_t Wsize = (argc > 4)? uint32_t(atoi(argv[4])) : 10;
	// OWM overlap
	real_t overlap = (argc > 5)? atof(argv[5]) : 0.8;

	// bintree_test(inputTXT, maxNumber);
	// bintree_traverse(inputTXT, maxNumber);
	// bintree64_test(inputTXT, maxNumber);
	// bintree64_traverse(inputTXT, maxNumber);
	// ordered_test(inputTXT, maxNumber);
	// ordered_traverse(inputTXT, maxNumber);
	// octree_test(inputTXT, maxNumber);
#ifndef HETER
	octree_traverse(inputTXT, maxNumber, Wsize, overlap);
#else
	octree_traverse_heter(inputTXT, maxNumber, factor, Wsize, overlap);
#endif
	// sycl_test(inputTXT, maxNumber);
	// sycl_traverse(inputTXT, maxNumber);	
	return 0;
}
