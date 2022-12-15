#ifndef OCTREE_BUILDER_H

#define OCTREE_BUILDER_H

#pragma once

#include <assert.h>

// #include <bintree/cuda/bintree_builder.h>
#include <sycl_cuda_ordered/source/sycl_builder.h>
#include <octree/octree.h>
#include <basic/utils.h>
#include <basic/morton.h>

#include <thrust/device_vector.h>

struct Octree_builder {

    // Bintree_builder bintree;
    SYCL_builder bintree;

    cl::sycl::queue& device_queue;

    // thrust::device_vector<aabb_t>           m_aabb;
    // thrust::device_vector<Split_task>       m_task_queues[2];
    // thrust::device_vector<uint32_t>         m_counters;
    // thrust::device_vector<octree_node>      m_octree;
    
    sycl::buffer<uint32_t, 1> b_counters;
    sycl::buffer<octree_node, 1> b_octree;
    sycl::buffer<aabb_t, 1> b_aabb;

    sycl::buffer<Split_task, 1> b_task_queue1;
    sycl::buffer<Split_task, 1> b_task_queue2;

    // Split_task* m_task_queues[2];
    
    // aabb_t*           m_aabb;
    // Split_task*       m_task_queues[2];
    // uint32_t*          m_counters;
    // octree_node*      m_octree;
    
    // uint32_t m_levels[64];
    uint32_t m_node_count;
    uint32_t m_leaf_count;
    uint32_t maxNumTasks;
    uint32_t maxNumNodes;

    Octree_builder( std::string inputTXT, const uint32_t chunk, cl::sycl::queue& q,
                    const uint32_t no, const uint32_t nn, const uint32_t nin, 
                    const uint32_t nt, const uint32_t non);

    void build();

    void reset();

    ~Octree_builder() {
        // free(m_aabb, device_queue);
        // free(m_octree, device_queue);
        // free(m_counters, device_queue);
        // free(m_task_queues[0], device_queue);
        // free(m_task_queues[1], device_queue);
    }
};

#include <sycl_cuda_octree/cuda/octree_builder.inl>

#endif // OCTREE_BUILDER_H