#ifndef OCTREE_BUILDER_H

#define OCTREE_BUILDER_H

#pragma once

#include <assert.h>

// #include <bintree/cuda/bintree_builder.h>
#include <ordered/cuda/ordered_builder.h>
#include <octree/octree.h>
#include <basic/utils.h>
#include <basic/morton.h>

#include <thrust/device_vector.h>

struct Octree_builder {

    // Bintree_builder bintree;
    Ordered_builder bintree;

    thrust::device_vector<aabb_t>           m_aabb;
    thrust::device_vector<Split_task>       m_task_queues[2];
    thrust::device_vector<uint32_t>         m_counters;
    thrust::device_vector<octree_node>      m_octree;

    aabb_t* b_aabb;
    // Split_task* b_task_queues[2];
    // uint32_t* b_counters;
    octree_node* b_octree;
    
    // uint32_t m_levels[64];
    uint32_t m_node_count;
    uint32_t m_leaf_count;
    uint32_t maxNumTasks;
    uint32_t maxNumNodes;

    Octree_builder( std::string inputTXT, const uint32_t chunk);

    void build();

    void reset();

    ~Octree_builder() {
        cudaFree(b_aabb);
        cudaFree(b_octree);
        // cudaFree(b_task_queues[0]);
        // cudaFree(b_task_queues[1]);
        // cudaFree(b_counters);
    }
};

#include <octree/cuda/octree_builder.inl>

#endif // OCTREE_BUILDER_H