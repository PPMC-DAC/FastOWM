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

    static const uint32_t kInvalid = uint32_t(-1);

    /// empty constructor
    inline
    octree_node() {}

    /// leaf constructor
    inline 
    octree_node(const uint32_t n, const uint32_t min, const uint32_t obj) : 
        bitmask(0u), first_child_id(kInvalid), numPts(n), min_idx(min), object_idx(obj) {} 

    /// full constructor
    inline 
    octree_node(const uint32_t mask, const uint32_t index, const uint32_t n, const uint32_t min) : 
        bitmask(mask), first_child_id( index ), numPts(n), min_idx(min), object_idx(kInvalid) {}

    // get the index of the i-th child (among the active ones)
    inline
    uint32_t get_child(const uint32_t i) const
    {
        if( bitmask & (1u << i) )
            return first_child_id + i;
        else
            return kInvalid;

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