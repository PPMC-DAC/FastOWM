//============================================================================
// Name			: ValueIterationScheduler.h
// Author		: Denisa Constantinescu
// Version		: 1.0
// Date			:
// Copyright	: Department. Computer's Architecture (c)
// Description	: ValueIterationScheduler.h implementation for ValueIteration project
//                The workload will be distributed row-wise, not by T matrix cell,
//                ({begin, end} in 0..NS)
//============================================================================

#ifndef LIDARSCHEDULERONEAPI_H
#define LIDARSCHEDULERONEAPI_H

//#define POLICY_IMPROVMENT
//#define DYNAMIC
//#define ORACLE

// #include "vi.h"

#include <cstdlib>
#include <iostream>

// #include <tbb/tick_count.h>
// #include <tbb/blocked_range.h>
// #include "tbb/parallel_for.h"
// #include "tbb/task.h"

// #include "../NavigationMDP/Common.h"

#include "DynamicOneApi.h"

#include "CL/sycl.hpp"
#include "Functors.h"

using namespace cl::sycl;

using namespace std;
using namespace tbb;

/*****************************************************************************
 * Global variables
 * **************************************************************************/

extern cl::sycl::queue gpu_queue;
extern cl::sycl::context ctx;
extern cl::sycl::event e1, e2, e3, e4, e5;

// extern class vi vi;
extern float *probability;
extern ulong *nextCellPos;
extern ulong *nextStatePos;
extern float *R;
//extern float*   s_Q;
//extern float*   s_V;

static size_t WORK_GROUP_SIZE = 256;

class LidarSchedulerOneApi
{
public:
#ifdef L_DEBUG
    uint32_t gpu_count;
    uint32_t cpu_count;
#endif

public:
    LidarSchedulerOneApi(/*class vi* h_v*/)
    {
#ifdef L_DEBUG
        gpu_count = 0;
        cpu_count = 0;
#endif        
        cout << "Instanciated LidarSchedulerOneApi \n";
    }



    /*This function launches the kernel*/
    void OperatorGPU(ulong begin, ulong end, cl::sycl::event *event)
    {
#ifdef L_DEBUG
        gpu_count += end - begin;
#endif

        // std::cout << "\t\tTrabajo para GPU: " << end - begin << " celdas\n" << std::flush;

        *event = gpu_queue.submit([&](cl::sycl::handler &cgh) {

#ifdef INDEX
            LidarEvaluationFunctorIdx kernel(  begin,
                                            initX,
                                            initY,
                                            nCols,
                                            Wsize,
                                            Displace,
                                            array_indexes,
                                            array_all_points,
                                            minNumPoints,
                                            minIDs);
#else
            LidarEvaluationFunctor kernel(  begin,
                                            initX,
                                            initY,
                                            nCols,
                                            Wsize,
                                            Displace,
                                            aligned_qtree,
                                            minNumPoints,
                                            minIDs);
#endif

            cgh.parallel_for(cl::sycl::range<1>(end - begin), kernel);
        });
    }


    /*Serial version of the code */
    void OperatorCPU(ulong begin, ulong end)
    {

#ifdef L_DEBUG
        cpu_count += end - begin;
#endif

        // std::cout << "\tTrabajo para CPU: " << end - begin << " celdas\n" << std::flush;

        // for(ulong i = begin; i<end; i++ ) {

        //   int cellPoints;

        //   Vector2D cellCenter = {initX + (i%nCols)*Displace, initY + (int)(i/nCols)*Displace};

        //   Lpoint newmin = gpuSearchNeighborsMin(cellCenter, aligned_qtree, Wsize*0.5, Wsize*0.5, cellPoints);

        //   if(cellPoints >= minNumPoints ){

        //     minIDs[i] = newmin.id;

        //   }

        // }

        Lpoint newmin = {0,0.0,0.0,std::numeric_limits<double>::max()};
        Vector2D cellCenter;
        Vector2D boxMax, boxMin;
        int cellPoints = 0;

        for(int idx = begin ; idx < end ; idx++ ){

            Vector2D cellCenter = {initX + (idx%nCols)*Displace, initY + (int)(idx/nCols)*Displace};

            makeBox(cellCenter, Wsize*0.5, boxMin, boxMax);

            if(insideBox2D(&newmin,boxMin,boxMax)){

                Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};

                int old_cellPoints = cellPoints;

#ifdef INDEX
                Lpoint tmp = gpuSearchIndex(oCell, array_indexes, array_all_points, Displace*0.5, Wsize*0.5, cellPoints);
#else                
                Lpoint tmp = gpuSearchNeighborsMin(oCell, aligned_qtree, Displace*0.5, Wsize*0.5, cellPoints);
#endif                

                // We're assuming the points were equidistant throughout the cell, which isn't always true.
                cellPoints += static_cast<int>(old_cellPoints * Overlap);

                if(tmp.z - newmin.z < 1e-5){
                    newmin = tmp;
                }

            } else {

#ifdef INDEX
                newmin = gpuSearchIndex(cellCenter, array_indexes, array_all_points, Wsize*0.5, Wsize*0.5, cellPoints);
#else                
                newmin = gpuSearchNeighborsMin(cellCenter, aligned_qtree, Wsize*0.5, Wsize*0.5, cellPoints);
#endif

            }

            if( cellPoints >= minNumPoints ){

                minIDs[idx] = newmin.id;

            }

        }
    }
};

#endif //LIDARSCHEDULERONEAPI_H
