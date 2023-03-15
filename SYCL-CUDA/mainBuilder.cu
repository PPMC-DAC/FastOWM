// #include <bintree/bintree_test.h>
// #include <bintree64/bintree64_test.h>
//#include <ordered/ordered_test.h>
#include <octree/octree_test.h>

int main(int argc, char* argv[])
{
	std::string inputTXT = (argc > 1)? argv[1] : "data/INAER_2011_Alcoy.xyz";
	uint32_t maxNumber = (argc > 2)? uint32_t(atoi(argv[2])) : 16;
	float factor = (argc > 3)? atof(argv[3]) : 0.75;

	// bintree_test(inputTXT, maxNumber);
	// bintree_traverse(inputTXT, maxNumber);
	// bintree64_test(inputTXT, maxNumber);
	// bintree64_traverse(inputTXT, maxNumber);
	// ordered_test(inputTXT, maxNumber);
	// ordered_traverse(inputTXT, maxNumber);
	// ordered_traverse_S3CPU(inputTXT, maxNumber);
	// octree_test(inputTXT, maxNumber);
	 octree_traverse(inputTXT, maxNumber);
	// octree_traverse_heter(inputTXT, maxNumber, factor);
	// sycl_test(inputTXT, maxNumber);
	// sycl_traverse(inputTXT, maxNumber);	
	return 0;
}
