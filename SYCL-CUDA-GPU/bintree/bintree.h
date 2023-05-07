#ifndef BINTREE_H

#define BINTREE_H

#pragma once

#include <basic/types.h>

struct bintree_node
{
    uint32_t parent_idx; // parent node
    uint32_t left_idx;   // index of left  child node
    uint32_t right_idx;  // index of right child node
    uint32_t numPts;     // number of points
    uint32_t min_idx;    // index of subtree's min point
    uint32_t object_idx; // == 0xFFFFFFFFu if internal node.
};

bintree_node default_node ={ 0xFFFFFFFFu,   // parent
                             0xFFFFFFFFu,   // left
                             0xFFFFFFFFu,   // right
                             0xFFFFFFFFu,   // numPts
                             0xFFFFFFFFu,   // min
                             0xFFFFFFFFu }; // object

#endif // BINTREE_H