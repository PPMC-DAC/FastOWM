#include <iostream>

#include <bintree/cuda/bintree_builder.h>

void bintree_test(std::string inputTXT, const uint32_t chunkDim)
{

	Bintree_builder builder(inputTXT, chunkDim);

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
      
    cudaEvent_t start, stop;
    cudaEventCreate( &start );
    cudaEventCreate( &stop );

    double time = 0.0;
#ifndef DEBUG
    int n_tests = 10;
#else
    int n_tests = 1;
#endif

    for (int i = 0; i < n_tests; i++)
    {
        builder.prefetchAll();
        
        std::chrono::time_point<tempo_t> start = tempo_t::now();
      
        builder.build();
    
        double dtime = cast_t(tempo_t::now() - start).count();

        time += dtime;
        std::cout << "Test: " << dtime << " ms\n";

        builder.reset();
    }
	time /= (float)(n_tests);
    
    
    std::cout << "  CREATION takes: " << time << " ms\n";

    // builder.prefetchAll();

    // builder.build();

    // point_t p = builder.point_cloud[0];
    // printf("%lf %lf %lf\n", p.x,p.y,p.z);
    // p.x -= builder.BBox.lower.x;
    // p.y -= builder.BBox.lower.y;
    // p.x /= builder.diffBox.x;
    // p.y /= builder.diffBox.y;
    // printf("%lf %lf %lf\n", p.x,p.y,p.z);
    // std::cout << uint32_t(p.x*65536.0) << "\n";
    // std::cout << uint32_t(p.y*65536.0) << "\n";
    // uint32_t mortonP = bintree_builder::morton2D(p);
    // std::cout << std::bitset<32>(mortonP) << "\n";
    // uint32_t dux = bintree_builder::CompactBy1(mortonP >> 1);
    // uint32_t duy = bintree_builder::CompactBy1(mortonP);
    // std::cout << std::bitset<32>(dux) << "\n";
    // std::cout << std::bitset<32>(duy) << "\n";
    // std::cout << dux << "\n";
    // std::cout << duy << "\n";
    // double dx = double(dux) / 65536.0;
    // double dy = double(duy) / 65536.0;
    // printf("%lf %lf %lf\n", dx,dy,p.z);
    // dx *= builder.diffBox.x;
    // dy *= builder.diffBox.y;
    // dx += builder.BBox.lower.x;
    // dy += builder.BBox.lower.y;
    // printf("%lf %lf %lf\n", dx,dy,p.z);

	return;
}

int main(int argc, char* argv[])
{

	std::string inputTXT = (argc > 1)? argv[1] : "../data/INAER_2011_Alcoy.xyz";
    uint32_t chunkDim = (argc > 2)? uint32_t(atoi(argv[2])) : 16;

	bintree_test(inputTXT, chunkDim);

	return 0;
}