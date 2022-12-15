#ifndef SYCL_BUILDER_H

#define SYCL_BUILDER_H

#pragma once

// #include <CL/sycl.hpp>
// #include <CL/sycl/backend/cuda.hpp>

#include <assert.h>

#include <basic/types.h>
#include <bintree/bintree.h>
#include <basic/utils.h>
#include <basic/morton.h>

using node_t = bintree_node;

struct SYCL_builder {

    uint32_t leafSize;
    
    uint32_t numObjects;
    uint32_t numLeafs;
    uint32_t numInternalNodes;
    uint32_t numNodes;

    whole_t BBox; // whole scenario
    point_t diffBox;

    cl::sycl::queue& device_queue;

    // int deviceId;
    int numberOfBlocks;
    int threadsPerBlock;

	sycl::buffer<point_t, 1>      point_cloud;
	sycl::buffer<point_t, 1>      ord_point_cloud;
	sycl::buffer<uint32_t, 1>      morton;
	sycl::buffer<uint32_t, 1>      indices;
	sycl::buffer<node_t, 1>      node_list;
	sycl::buffer<aabb_t, 1>      aabb_list;
    sycl::buffer<uint32_t, 1>      flags;

	// point_t*      point_cloud;
	// point_t*      ord_point_cloud;
	// uint32_t*      morton;
	// uint32_t*      indices;
	// node_t*      node_list;
	// aabb_t*     aabb_list;
    // uint32_t*      flags;

    SYCL_builder( std::string inputTXT, const uint32_t chunk, cl::sycl::queue& q,
                    const uint32_t no, const uint32_t nn, const uint32_t nin);

    // SYCL_builder( std::string inputTXT, const uint32_t chunk, cl::sycl::queue q);

    void build();

    void reset();
    
    /*destructor vacío porque al utilizar buffers no necesito liberar explicitamente memoria*/
    ~SYCL_builder(){}

    // ~SYCL_builder()
    // {

    // free(point_cloud, device_queue);
    // free(ord_point_cloud, device_queue);
    // free(morton, device_queue);
    // free(indices, device_queue);
    // free(node_list, device_queue);
    // free(aabb_list, device_queue);
    // free(flags, device_queue);

    // } // destructor
};

#include <sycl_cuda_ordered/source/sycl_builder.inl>

#endif // SYCL_BUILDER_H