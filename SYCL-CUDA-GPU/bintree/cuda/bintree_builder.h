#ifndef BINTREE_BUILDER_H

#define BINTREE_BUILDER_H

#pragma once

#include <assert.h>

#include <bintree/bintree.h>
#include <basic/utils.h>
#include <basic/morton.h>

struct Bintree_builder {

	point_t*        point_cloud;
	uint32_t*       morton;
	uint32_t*       morton_out;
	uint32_t*       indices;
	uint32_t*       indices_out;
	// uint64_t*       morton64;
	bintree_node*   node_list;
	aabb_t*         aabb_list;
    uint32_t*       flags;
    void*           d_temp_storage;

    uint32_t leafSize;
    
    uint32_t numObjects;
    uint32_t numLeafs;
    uint32_t numInternalNodes;
    uint32_t numNodes;

    whole_t BBox; // whole scenario
    point_t diffBox;

	int deviceId;
	int numberOfSMs;
    int numberOfBlocks;
    int threadsPerBlock;

    Bintree_builder( std::string inputTXT, const uint32_t chunk);

    void build();

    void prefetch();

    template <typename pointer_t, typename size_t, typename ...Rest>
    void prefetch( pointer_t p, size_t s, int device, Rest... args );

    void prefetchAll();

    void reset();

    template<typename kernel_t>
    int getOccupancy( kernel_t kernel );

    ~Bintree_builder()
    {
        cudaFree(point_cloud);
        cudaFree(morton);
        cudaFree(morton_out);
        cudaFree(indices);
        cudaFree(indices_out);
        // cudaFree(morton64);
        cudaFree(node_list);
        cudaFree(aabb_list);
        cudaFree(flags);
        cudaFree(d_temp_storage);
    }
};

#include <bintree/cuda/bintree_builder.inl>

#endif // BINTREE_BUILDER_H