#ifndef OCTREE_BUILDER_H

#define OCTREE_BUILDER_H

#pragma once

#include <assert.h>

// #include <bintree/cuda/bintree_builder.h>
#include <sycl_ordered/source/sycl_builder.h>
#include <sycl_octree/octree.h>
// #include <basic/utils_sycl.h>
// #include <basic/morton_sycl.h>

// #include <thrust/device_vector.h>

struct Octree_builder {

    // Bintree_builder bintree;
    SYCL_builder bintree;

    // thrust::device_vector<aabb_t>           m_aabb;
    // thrust::device_vector<Split_task>       m_task_queues[2];
    // thrust::device_vector<uint32_t>         m_counters;
    // thrust::device_vector<octree_node>      m_octree;

    sycl::queue& device_queue;

    aabb_t* m_aabb;
    Split_task* m_task_queues[2];
    uint32_t* m_counters;
    octree_node* m_octree;
    
    // uint32_t m_levels[64];
    uint32_t m_node_count;
    uint32_t m_leaf_count;
    uint32_t maxNumTasks;
    uint32_t maxNumNodes;

    Octree_builder( std::string inputTXT, const uint32_t chunk, sycl::queue& q);

    void build();

    void reset();

    ~Octree_builder() 
    {
        free(m_aabb, device_queue);
        free(m_task_queues[0], device_queue);
        free(m_task_queues[1], device_queue);
        free(m_counters, device_queue);
        free(m_octree, device_queue);
    }
};

#include <sycl_octree/source/octree_builder.inl>

#endif // OCTREE_BUILDER_H