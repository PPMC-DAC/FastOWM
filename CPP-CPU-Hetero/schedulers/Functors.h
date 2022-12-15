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
    double initX;
    double initY;
    uint32_t nCols;
    uint32_t Wsize;
    double Displace;
    QtreeG4 qtree;
    uint32_t minNumPoints;
    int* minIDs;

public:
    LidarEvaluationFunctor( uint32_t offset_,
                            double initX_,
                            double initY_,
                            uint32_t nCols_,
                            uint32_t Wsize_,
                            double Displace_,
                            QtreeG4 qtree_,
                            uint32_t minNumPoints_,
                            int* minIDs_) : offset(offset_), initX(initX_), initY(initY_), nCols(nCols_), Wsize(Wsize_),
                                            Displace(Displace_), qtree(qtree_), minNumPoints(minNumPoints_),
                                            minIDs(minIDs_) {}

    void operator()(cl::sycl::id<1> index) const
    {
        uint32_t i = static_cast<int>(index[0]) + offset;

        int cellPoints;

        Vector2D cellCenter = {initX + (i%nCols)*Displace, initY + (int)(i/nCols)*Displace};

        Lpoint newmin = gpuSearchNeighborsMin(cellCenter, qtree, Wsize*0.5, Wsize*0.5, cellPoints);

        if(cellPoints >= minNumPoints ){
            minIDs[i] = newmin.id;
        }
    }

};



class LidarEvaluationFunctorIdx
{
private:
    uint32_t offset;
    double initX;
    double initY;
    uint32_t nCols;
    uint32_t Wsize;
    double Displace;
    QtreeG5 qtree;
    Lpoint* points;
    uint32_t minNumPoints;
    int* minIDs;

public:
    LidarEvaluationFunctorIdx( uint32_t offset_,
                            double initX_,
                            double initY_,
                            uint32_t nCols_,
                            uint32_t Wsize_,
                            double Displace_,
                            QtreeG5 qtree_,
                            Lpoint* points_,
                            uint32_t minNumPoints_,
                            int* minIDs_) : offset(offset_), initX(initX_), initY(initY_), nCols(nCols_), Wsize(Wsize_),
                                            Displace(Displace_), qtree(qtree_), points(points_), minNumPoints(minNumPoints_),
                                            minIDs(minIDs_) {}

    void operator()(cl::sycl::id<1> index) const
    {
        uint32_t i = static_cast<int>(index[0]) + offset;

        int cellPoints;

        Vector2D cellCenter = {initX + (i%nCols)*Displace, initY + (int)(i/nCols)*Displace};

        Lpoint newmin = gpuSearchIndex(cellCenter, qtree, points, Wsize*0.5, Wsize*0.5, cellPoints);

        if(cellPoints >= minNumPoints ){
            minIDs[i] = newmin.id;
        }
    }
};


#endif //FUNCTORS_H
