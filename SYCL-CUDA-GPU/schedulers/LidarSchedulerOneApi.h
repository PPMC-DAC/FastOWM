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

#include "sycl/sycl.hpp"
#include "Functors.h"

using namespace sycl;

using namespace std;
using namespace tbb;

/*****************************************************************************
 * Global variables
 * **************************************************************************/

extern Octree_builder* g_builder;
extern LBVHoct* g_lbvh_cpu;

extern sycl::queue device_queue;
// extern sycl::context ctx;
// extern sycl::event e1, e2, e3, e4, e5;

// extern class vi vi;
// extern float *probability;
// extern ulong *nextCellPos;
// extern ulong *nextStatePos;
// extern float *R;
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
    LBVHoct lbvh_cpu;

public:
    LidarSchedulerOneApi(/*class vi* h_v*/ const LBVHoct& lbvh) : lbvh_cpu(lbvh)
    {
#ifdef L_DEBUG
        gpu_count = 0;
        cpu_count = 0;
#endif        
        cout << "Instanciated LidarSchedulerOneApi \n";
    }



    /*This function launches the kernel*/
    void OperatorGPU(ulong begin, ulong end, sycl::event *event)
    {
#ifdef L_DEBUG
        gpu_count += end - begin;
#endif

        // const size_t total_size = (end - begin) * sizeof(uint32_t);
        // const size_t chunk_size = 1024 * 1024; // 1MB chunks, adjust as needed

        // for (size_t offset = 0; offset < total_size; offset += chunk_size) {
        //     size_t size = std::min(chunk_size, total_size - offset);
        //     device_queue.prefetch(minIDs + offset / sizeof(uint32_t), size);
        // }
        
        // Prefetch the exact chunk of minIDs we'll be working with
        device_queue.prefetch(minIDs + begin, (end - begin) * sizeof(uint32_t));

        *event = device_queue.submit([&](sycl::handler &cgh) {

            LidarEvaluationFunctor kernel(  begin,
                                            initBox,
                                            nRows,
                                            nCols,
                                            // Wsize,
                                            Displace,
                                            g_builder,
                                            minNumPoints,
                                            minIDs );

            cgh.parallel_for(
                sycl::nd_range<1>(sycl::range<1>(end - begin), sycl::range<1>(WORK_GROUP_SIZE)),
                kernel
            );

        });
    }


    /*Serial version of the code */
    void OperatorCPU(ulong begin, ulong end)
    {

// #ifdef L_DEBUG
//         cpu_count += end - begin;
// #endif

        aabb_t cellBox, oCell;
        uint32_t pointsCount, old_cellPoints, tmpId;
        uint32_t idmin = 0xFFFFFFFFu;

        for(int idx = begin ; idx < end ; idx++ ){
            
            int ii = idx % nCols;
            int jj = idx / nCols;

            cellBox.upper.y = jj*Displace + initBox.upper.y;
            cellBox.lower.y = jj*Displace + initBox.lower.y;
            cellBox.upper.x = ii*Displace + initBox.upper.x;
            cellBox.lower.x = ii*Displace + initBox.lower.x;

            if(idmin != 0xFFFFFFFFu && insideBox(cellBox, lbvh_cpu.point_list[idmin]))
            {

                // Vector2D oCell = {cellCenter.x + Wsize*0.5 - Displace*0.5 , cellCenter.y};
                oCell.upper.x = cellBox.upper.x;
                oCell.lower.x = cellBox.upper.x - Displace;
                oCell.upper.y = cellBox.upper.y;
                oCell.lower.y = cellBox.lower.y;

                old_cellPoints = pointsCount;

                pointsCount = 0u;
                tmpId = 0xFFFFFFFFu;

                traverseIterativeCPU(lbvh_cpu, oCell, pointsCount, tmpId);

                // We're assuming the points were equidistant throughout the cell, which isn't always true.

                /*En este punto, si queremos ser estrictos, en vez de hacer esta suposición podemos
                lanzar una búsqueda "countMin" en la zona en la que conocemos el mínimo, pudiendo lanzar
                las dos búsquedas diferentes en paralelo con un "parallel_invoke" */
                pointsCount += (uint32_t)(old_cellPoints * Overlap);

                /*Si he hecho traverseIterativeCPU y tmpId = 0xFFFFFFFFu, quiere decir que la BB con
                la que analizo no solapa ninguna BB de un nodo, lo que puede suceder muy facilmente
                porque la ventana de análisis tiene un margen de tamaño (Wsize-Displace) por cada
                lateral de la nube de puntos*/
                if(tmpId != 0xFFFFFFFFu && lbvh_cpu.point_list[tmpId].z < lbvh_cpu.point_list[idmin].z)
                {
                    idmin = tmpId;
                }

            } else {

                pointsCount = 0u;
                idmin = 0xFFFFFFFFu;

                traverseIterativeCPU(lbvh_cpu, cellBox, pointsCount, idmin);

            }

            minIDs[idx] = ( minNumPoints <= pointsCount )? idmin : 0u;

        }
    }
};

#endif //LIDARSCHEDULERONEAPI_H
