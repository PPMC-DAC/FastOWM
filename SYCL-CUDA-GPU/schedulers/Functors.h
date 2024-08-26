//
// Created by root on 4/1/20.
//

#ifndef FUNCTORS_H
#define FUNCTORS_H

#include "CL/sycl.hpp"

class LidarEvaluationFunctor
{
private:
    uint32_t offset;
    aabb_t initBox;
    LBVHoct lbvh;
    uint32_t nRows;
    uint32_t nCols;
    // uint32_t Wsize;
    double Displace;
    uint32_t minNumPoints;
    uint32_t* minIDs;

public:
    LidarEvaluationFunctor( uint32_t offset_,
                            aabb_t initBox_,
                            uint32_t nRows_,
                            uint32_t nCols_,
                            // uint32_t Wsize_,
                            double Displace_,
                            Octree_builder* builder_,
                            uint32_t minNumPoints_,
                            uint32_t* minIDs_) :    offset(offset_), 
                                                    nRows(nRows_),  
                                                    nCols(nCols_),
                                                    Displace(Displace_),
                                                    minNumPoints(minNumPoints_),
                                                    minIDs(minIDs_),
                                                    lbvh(builder_->m_octree, builder_->m_aabb, builder_->bintree.ord_point_cloud),
                                                    initBox(initBox_) {}

    void operator()(sycl::nd_item<1> item) const
    {
        const int id = item.get_global_id(0) + offset;

        const int idx = id % nCols;
        const int jdx = id / nCols;
        
        aabb_t cellBox;
        cellBox.upper.x = idx*Displace + initBox.upper.x;
        cellBox.lower.x = idx*Displace + initBox.lower.x;
        cellBox.upper.y = jdx*Displace + initBox.upper.y;
        cellBox.lower.y = jdx*Displace + initBox.lower.y;

        uint32_t pointsCount = 0;
        uint32_t idmin = 0xFFFFFFFFu;

        traverseIterative(lbvh, cellBox, pointsCount, idmin);

        minIDs[id] = ( minNumPoints <= pointsCount )? idmin : 0u;

        // if (minNumPoints <= pointsCount)
        // {
        //     minIDs[id] = idmin;
        // }

        // cosnt uint32_t result = (minNumPoints <= pointsCount) ? idmin : 0u;

        // // Asynchronous write to global memory
        // sycl::atomic_ref<uint32_t, sycl::memory_order::relaxed, 
        //                  sycl::memory_scope::device, 
        //                  sycl::access::address_space::global_space> 
        // atomic_minIDs(minIDs[id]);

        // atomic_minIDs.store(result);

        return;
    }

};

#endif //FUNCTORS_H
