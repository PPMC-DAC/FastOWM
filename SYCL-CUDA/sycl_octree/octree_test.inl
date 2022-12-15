#include <iostream>

#include <sycl_octree/source/octree_builder.h>

void octree_test(std::string inputTXT, const uint32_t chunkDim)
{

#ifdef FEMU
    cl::sycl::INTEL::fpga_emulator_selector device_selector{};
#elif GPU
    cl::sycl::gpu_selector device_selector{};
#elif CPU
    cl::sycl::cpu_selector device_selector{};
#endif

    cl::sycl::queue device_queue(device_selector);

    std::cout << "Device : " << device_queue.get_device().get_info<sycl::info::device::name>()  
    << " @ " << device_queue.get_device().get_info<sycl::info::device::max_clock_frequency>() << "Mhz (" << 
    device_queue.get_device().get_info<sycl::info::device::max_compute_units>() << " cores)" << std::endl;
    std::cout << device_queue.get_device().get_info<sycl::info::device::max_work_item_dimensions>() << " max_dim" << std::endl;
    // std::cout << device_queue.get_device().get_info<sycl::info::device::max_work_item_sizes>() << " max_size" << std::endl;
    std::cout << device_queue.get_device().get_info<sycl::info::device::max_work_group_size>() << " max_group" << std::endl;

	Octree_builder builder(inputTXT, chunkDim, device_queue);

#ifdef DEBUG
    std::cout << inputTXT << "; " << builder.bintree.numObjects << " points; chunkDim: " << builder.bintree.leafSize << std::endl;
    std::cout << builder.bintree.numLeafs << " leaf nodes\n";
    std::cout << builder.bintree.numInternalNodes << " internal nodes\n";
    std::cout << builder.bintree.numNodes << " nodes\n\n";

    std::cout << "BBox: " << builder.bintree.BBox.upper.x << "," << builder.bintree.BBox.upper.y << " ";
    std::cout << builder.bintree.BBox.lower.x << "," << builder.bintree.BBox.lower.y << "\n\n";
  
    std::cout << "diffBox: " << builder.bintree.diffBox.x << "," << builder.bintree.diffBox.y << "\n\n";
#endif

    double time = 0.0;
#ifndef DEBUG
    int n_tests = 100;
    std::cout << "Performing " << n_tests << " tests:\n";
#else
    int n_tests = 1;
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
    
    
    std::cout << "  OCTREE creation takes: " << time << " ms\n";

#ifdef DEBUG
    printf("  points     : %u\n", builder.bintree.numObjects);
    printf("  points/sec : %f M\n", (builder.bintree.numObjects / time) / 1.0e6f );
    printf("  bintree nodes : %u\n", builder.bintree.numNodes);
    printf("  octree nodes  : %u\n", builder.m_node_count );
    printf("  leaves : %u\n", builder.m_leaf_count );
    // for (uint32_t level = 0; level < 16; ++level)
    //     fprintf(stderr, "  level %u : %u nodes\n", level, builder.m_levels[level+1] - builder.m_levels[level] );
#endif

    // const octree_node root = builder.m_octree[3];
    // for(int i=0; i<8; i++){
    //     printf("Hijo %d: %u\n", i, root.get_child(i));
    // }

    // uint32_t hojas = 0;
    // uint32_t primera_hoja = 0;
    // for(int i=0; i<builder.m_node_count; i++)
    // {
    //     const octree_node node = builder.m_octree[i];
    //     // if(node.bitmask != 0u) {
    //     //     printf("[%d] :: ", i);
    //     //     for(int j=0; j<8; j++)
    //     //         printf("%u ", node.get_child(j) );
    //     //     printf("\n\n");
    //     // }
        
    //     if(node.bitmask == 0u) hojas++; 

    //     if(node.bitmask == 0u && primera_hoja == 0) primera_hoja=i;

    //     // uint32_t active_child = 0;
    //     // for(int j=0; j<8; j++)
    //     //     if(node.get_child(j) != uint32_t(-1)) active_child++;
    //     // if(active_child != 8 && active_child != 0)
    //     //     printf("Hay un nodo INCOMPLETO\n");
    // }
    // printf("Tengo %u hojas\n", hojas);

    // for(int i=primera_hoja; i<builder.m_node_count; i++)
    // {
    //     const octree_node node = builder.m_octree[i];
    //     if(node.bitmask != 0u){
    //         printf("Hay nodos MEZCLADOS!!\n");
    //         break;
    //     }
    // }

	return;
}
