#ifndef SYCL_BUILDER_H

#define SYCL_BUILDER_H

#pragma once

#include <assert.h>

#include <basic/types_sycl.h>
#include <bintree/bintree.h>
#include <basic/utils_sycl.h>
#include <basic/morton_sycl.h>

#include <sycl/sycl.hpp>

struct SYCL_builder {

	point_t*        point_cloud;
	point_t*        ord_point_cloud;
	uint32_t*       morton;
	uint32_t*       morton_out;
	uint32_t*       indices;
	uint32_t*       indices_out;
	// uint64_t*       morton64;
	bintree_node*   node_list;
	aabb_t*         aabb_list;
    uint32_t*       flags;
    // void*           d_temp_storage;

    uint32_t leafSize;
    
    uint32_t numObjects;
    uint32_t numLeafs;
    uint32_t numInternalNodes;
    uint32_t numNodes;

    whole_t BBox; // whole scenario
    point_t diffBox;

    sycl::queue& device_queue;

    size_t num_groups;
    size_t items_per_group;
    size_t grid_dim;

    SYCL_builder( std::string inputTXT, const uint32_t chunk, sycl::queue& q);

    void build();

    void reset();

    ~SYCL_builder();
};

#include <sycl_ordered/source/sycl_builder.inl>

#endif // SYCL_BUILDER_H