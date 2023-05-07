#include <iostream>

#include <sycl_cuda_octree/cuda/octree_builder.h>

class CUDASelector : public sycl::device_selector {
public:
  int operator()(const sycl::device &Device) const override {
    using namespace sycl::info;

    const std::string DriverVersion = Device.get_info<device::driver_version>();

    if (Device.is_gpu() && (DriverVersion.find("CUDA") != std::string::npos)) {
      std::cout << "CUDA device found: ";
      return 1;
    };
    return -1;
  }
};

void octree_test(std::string inputTXT, const uint32_t chunkDim)
{

    cl::sycl::queue device_queue(CUDASelector().select_device());

    std::cout << device_queue.get_device().get_info<sycl::info::device::name>()  
    << " @ " << device_queue.get_device().get_info<sycl::info::device::max_clock_frequency>() << "Mhz (" << 
    device_queue.get_device().get_info<sycl::info::device::max_compute_units>() << " cores)" << std::endl;

    uint32_t numObjects = 0;
    whole_t BBox;

#ifndef CHECK
    readHeader(inputTXT, BBox, numObjects);
#else
    inputTXT = "CHECK";
    chunkDim = 2;
    numObjects = 32;
#endif

    uint32_t numLeafs = (uint32_t)((numObjects-1)/chunkDim) + 1u;
    uint32_t numInternalNodes = numLeafs - 1u;
    uint32_t numNodes = 2*numLeafs - 1u;

    // Numero minimo de nodos en bloques de 8 que contienen a todas las hojas
    uint32_t aux = uint32_t((numLeafs-1)/8) + 1u;
    // pero necesito todos los nodos intermedios también
    int i=0;
    uint32_t maxNumTasks = 0;
    uint32_t pot = 0;
    for(; pot<aux; i++) {
        pot = uint32_t(pow(8,i));
        maxNumTasks += pot;
    }

    /* Procesar un nodo hoja no va a ser una tarea para nosotros, 
    pero sí que va a ser un nuevo nodo del octree*/

    uint32_t maxNumNodes = maxNumTasks + numLeafs;

#ifdef DEBUG
    printf("Necesito %d niveles en el arbol\n", i);
    printf("Voy a tener como máximo %u nodos\n", maxNumNodes);
#endif

    // SYCL_builder builder(inputTXT, chunkDim, device_queue, numObjects, numNodes, numInternalNodes);

	Octree_builder builder(inputTXT, chunkDim, device_queue, numObjects, numNodes, numInternalNodes, maxNumTasks, maxNumNodes);

#ifdef DEBUG
    std::cout << inputTXT << "; " << builder.bintree.numObjects << " points; chunkDim: " << builder.bintree.leafSize;
    std::cout << "; threadsPerBlock: " << builder.bintree.threadsPerBlock << "; numBlocks: " << builder.bintree.numberOfBlocks << std::endl;
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

        builder.bintree.reset();
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
