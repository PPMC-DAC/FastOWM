
#include <oneapi/dpl/execution>
#include <oneapi/dpl/algorithm>

#include <sycl_octree/octree_test.h>

int main(int argc, char* argv[])
{
	std::string inputTXT = (argc > 1)? argv[1] : "data/INAER_2011_Alcoy.xyz";
    uint32_t maxNumber = (argc > 2)? uint32_t(atoi(argv[2])) : 16; // Number of LiDAR points per leaf-node
	int chunkGPU = (argc > 3)? atoi(argv[3]) : 1024; // Only used for octree_traverse_scheduler 
	// OWM window size
	uint32_t Wsize = (argc > 4)? uint32_t(atoi(argv[4])) : 10;
	// OWM overlap
	real_t overlap = (argc > 5)? atof(argv[5]) : 0.8;

	octree_traverse_scheduler(inputTXT, maxNumber, chunkGPU, Wsize, overlap);

	return 0;
}
