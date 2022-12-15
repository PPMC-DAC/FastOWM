#include <iostream>

#include <bintree64/cuda/bintree64_builder.h>

void bintree64_test(std::string inputTXT, const uint32_t chunkDim)
{

	Bintree64_builder builder(inputTXT, chunkDim);

#ifdef DEBUG
    std::cout << inputTXT << "; " << builder.numObjects << " points; chunkDim: " << builder.leafSize;
    std::cout << "; threadsPerBlock: " << builder.threadsPerBlock << "; numBlocks: " << builder.numberOfBlocks << std::endl;
    std::cout << builder.numLeafs << " leaf nodes\n";
    std::cout << builder.numInternalNodes << " internal nodes\n";
    std::cout << builder.numNodes << " nodes\n\n";

    std::cout << "BBox: " << builder.BBox.upper.x << "," << builder.BBox.upper.y << " ";
    std::cout << builder.BBox.lower.x << "," << builder.BBox.lower.y << "\n\n";
  
    std::cout << "diffBox: " << builder.diffBox.x << "," << builder.diffBox.y << "\n\n";
#endif

    double time = 0.0;
#ifndef DEBUG
    int n_tests = 100;
    std::cout << "Performing " << n_tests << " tests:\n";
#else
    int n_tests = 4;
    std::cout << "Performing " << n_tests << " test:\n";
#endif

    for (int i = 0; i <= n_tests; i++)
    {
        std::chrono::time_point<tempo_t> start = tempo_t::now();
      
        builder.build();
    
        double dtime = cast_t(tempo_t::now() - start).count();
        if(i)
            time += dtime;
#ifdef DEBUG
        std::cout << "Test: " << dtime << " ms\n";
#endif

        builder.reset();
    }
	time /= (float)(n_tests);
    
    
    std::cout << "  CREATION takes: " << time << " ms\n";


	return;
}
