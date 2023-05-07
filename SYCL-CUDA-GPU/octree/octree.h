#ifndef OCTREE_H

#define OCTREE_H

#pragma once

struct octree_node
{
    // uint32_t parent_idx; // parent node
    // uint32_t children[8];
    // uint32_t numPts;     // number of points
    // uint32_t min_idx;    // index of subtree's min point
    // uint32_t object_idx; // == 0xFFFFFFFFu if internal node.

    // static const uint32_t kInvalid = uint32_t(-1);

    /// empty constructor
    // __host__ __device__ __forceinline__
    // octree_node() {}

    /// leaf constructor
    __host__ __device__ __forceinline__ constexpr
    octree_node(const uint32_t n, const uint32_t min, const uint32_t obj) : 
        bitmask(0u), first_child_id(uint32_t(-1)), numPts(n), min_idx(min), object_idx(obj) {} 

    /// full constructor
    __host__ __device__ __forceinline__ constexpr 
    octree_node(const uint32_t mask, const uint32_t index, const uint32_t n, const uint32_t min) : 
        bitmask(mask), first_child_id( index ), numPts(n), min_idx(min), object_idx(uint32_t(-1)) {}

    // get the index of the i-th child (among the active ones)
    __host__ __device__ __forceinline__ constexpr
    uint32_t get_child(const uint32_t i) const
    {
        if( bitmask & (1u << i) )
            return first_child_id + i;
        else
            return uint32_t(-1);

        // return (bitmask & (1u << i))? first_child_id + i : kInvalid;
    }

    uint32_t bitmask;
    uint32_t first_child_id;
    // uint32_t aabb_idx;
    uint32_t numPts;
    uint32_t min_idx;
    uint32_t object_idx;
};


#endif // OCTREE_H