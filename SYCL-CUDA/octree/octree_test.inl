#include <iostream>

#include <octree/cuda/octree_builder.h>

void octree_test(std::string inputTXT, const uint32_t chunkDim)
{

	Octree_builder builder(inputTXT, chunkDim);

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

    for (int i = 0; i < n_tests; i++)
    {
        std::chrono::time_point<tempo_t> start = tempo_t::now();
      
        builder.build();
    
        double dtime = cast_t(tempo_t::now() - start).count();
        
        // if(i)
            time += dtime;
#ifdef DEBUG
        std::cout << "Test: " << dtime << " ms\n";
#endif

#ifdef DEBUG
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "parent_idx, left_idx, right_idx, object_idx, numPts, min_idx (Z)" << std::endl;
        std::cout << std::endl;
        for(int i = 0; i<builder.bintree.numNodes; i++){
            std::cout << builder.bintree.node_list[i].parent_idx << ", " << builder.bintree.node_list[i].left_idx << ", ";
            std::cout << builder.bintree.node_list[i].right_idx << ", " << builder.bintree.node_list[i].object_idx << ", ";
            std::cout << builder.bintree.node_list[i].numPts << ", " << builder.bintree.node_list[i].min_idx << "(";
            if(builder.bintree.node_list[i].min_idx != 0xFFFFFFFFu)
                std::cout << builder.bintree.ord_point_cloud[builder.bintree.node_list[i].min_idx].z;
            std::cout << ")\n";
        }
        std::cout << std::endl;
        std::cout << "upper.x, upper.y, lower.x, lower.y" << std::endl;
        std::cout << std::endl;
        for(int i = 0; i<builder.bintree.numNodes; i++){
            std::cout << builder.bintree.aabb_list[i].upper.x << "," << builder.bintree.aabb_list[i].upper.y << " , ";
            std::cout << builder.bintree.aabb_list[i].lower.x << "," << builder.bintree.aabb_list[i].lower.y << "\n";
        }
        std::cout << std::endl;
        std::cout << std::endl;
        std::cout << "bitmask, firs_child_id, object_idx, numPts, min_idx (Z)" << std::endl;
        std::cout << std::endl;
        for(int i = 0; i<builder.m_node_count+2; i++){
            std::cout << std::bitset<8>(builder.b_octree[i].bitmask) << ", ";
            std::cout << builder.b_octree[i].first_child_id << ", " << builder.b_octree[i].object_idx << ", ";
            std::cout << builder.b_octree[i].numPts << ", " << builder.b_octree[i].min_idx << "(";
            if(builder.b_octree[i].min_idx != 0xFFFFFFFFu)
                std::cout << builder.bintree.ord_point_cloud[builder.b_octree[i].min_idx].z;
            std::cout << ")\n";

            // std::cout << std::endl;
            // std::cout << std::endl;
            // std::cout << std::endl;
            std::cout << "Hijos de " << i << ": ";
            std::cout << builder.b_octree[i].get_child(0) << " ";
            std::cout << builder.b_octree[i].get_child(1) << " ";
            std::cout << builder.b_octree[i].get_child(2) << " ";
            std::cout << builder.b_octree[i].get_child(3) << " ";
            std::cout << builder.b_octree[i].get_child(4) << " ";
            std::cout << builder.b_octree[i].get_child(5) << " ";
            std::cout << builder.b_octree[i].get_child(6) << " ";
            std::cout << builder.b_octree[i].get_child(7) << " ";
            std::cout << std::endl;
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << "upper.x, upper.y, lower.x, lower.y" << std::endl;
        std::cout << std::endl;
        for(int i = 0; i<builder.m_node_count; i++){
            std::cout << builder.b_aabb[i].upper.x << "," << builder.b_aabb[i].upper.y << " , ";
            std::cout << builder.b_aabb[i].lower.x << "," << builder.b_aabb[i].lower.y << "\n";
        }
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
