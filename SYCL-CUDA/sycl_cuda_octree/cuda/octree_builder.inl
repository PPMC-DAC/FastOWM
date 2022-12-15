#include <algorithm>
// #include <cub/cub.cuh>
// #include <cub/device/device_radix_sort.cuh>

template <typename vector_type>
void need_space(vector_type& vec, const uint32_t size)
{
  if (vec.size() < size)
    vec.resize( size );
}

/// intra-warp inclusive scan
///
/// \param val      per-threrad input value
/// \param tidx     warp thread index
/// \param red      scan result storage (2*WARP_SIZE elements)
template <typename T> 
__device__ __forceinline__ T scan_warp(T val, const int32_t tidx, volatile T *red)
{
    // pad initial segment with zeros
    red[tidx] = 0;
    red += 32;

    // Hillis-Steele scan
    red[tidx] = val;
    val += red[tidx-1];  red[tidx] = val;
    val += red[tidx-2];  red[tidx] = val;
    val += red[tidx-4];  red[tidx] = val;
    val += red[tidx-8];  red[tidx] = val;
    val += red[tidx-16]; red[tidx] = val;
	return val;
}
/// return the total from a scan_warp
///
/// \param red      scan result storage
template <typename T> 
__device__ __forceinline__ T scan_warp_total(volatile T *red) { return red[63]; }

/// alloc n elements per thread from a common pool, using a synchronous warp scan
///
/// \param n                number of elements to alloc
/// \param warp_tid         warp thread index
/// \param warp_red         temporary warp scan storage (2*WARP_SIZE elements)
/// \param warp_broadcast   temporary warp broadcasting storage
__device__ __forceinline__
uint32_t myalloc(uint32_t n, uint32_t* pool, const int32_t warp_tid, volatile uint32_t* warp_red, volatile uint32_t* warp_broadcast)
{
    uint32_t warp_scan  = scan_warp( n, warp_tid, warp_red ) - n;
    uint32_t warp_count = scan_warp_total( warp_red );
    if (warp_tid == 0)
        *warp_broadcast = atomicAdd( pool, warp_count );

    return *warp_broadcast + warp_scan;
}



__global__
void init_tasks(Split_task* queue1, Split_task* queue2, const uint32_t num_tasks)
{
    uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
    uint32_t stride = blockDim.x * gridDim.x;

    for(uint32_t i = index; i<num_tasks; i += stride)
    {
        queue1[i] = Split_task(0u,0u,0u,0u);
        queue2[i] = Split_task(0u,0u,0u,0u);
    }

    return;
}

__global__
void init_octree(octree_node* nodes, aabb_t* bb, const uint32_t num_nodes, const aabb_t defaul_aabb)
{
    uint32_t index = threadIdx.x + blockIdx.x * blockDim.x;
    uint32_t stride = blockDim.x * gridDim.x;

    for(uint32_t i = index; i<num_nodes; i += stride)
    {
        nodes[i] = octree_node(0u,0u,0u,0u);
        bb[i] = defaul_aabb;
    }

    return;
}




namespace octree_builder {


/// Utility class to hold results of an octree collection step
struct Octree_collection
{
    uint32_t node_count;
    uint32_t leaf_count;
    uint32_t bitmask;
};

// Helper class to collect the 8 children of a given node from
// a binary kd-tree structure.
// The class implements template based compile-time recursion.
template <uint32_t LEVEL, uint32_t STRIDE>
struct Octree_collector
{
    static __device__ void find_children(
        const uint32_t        node_index,
        const bintree_node*      nodes,
        Octree_collection*  result,
        uint32_t*             children,
        const uint32_t        octant = 0)
    {
        const bintree_node node = nodes[ node_index ];

        const bool active0 = (node.left_idx == 0xFFFFFFFFu)? 0 : 1;
        const bool active1 = (node.right_idx == 0xFFFFFFFFu)? 0 : 1;

        if ((active0 == false) && (active1 == false))
        {
            // here we have some trouble... this guy is likely not
            // fitting within an octant of his parent.
            // Which octant do we assign this guy to?
            // Let's say the minimum...
            children[ STRIDE * result->node_count ] = node_index;
            result->bitmask |= 1u << (octant << (3 - LEVEL));
            result->node_count += 1;
            result->leaf_count += 1;
        }
        else
        {
            // traverse the children in Morton order: preserving
            // the order is key here, as we want the output counter
            // to match the bitmask pop-counts.
            if (active0)
            {
                Octree_collector<LEVEL+1,STRIDE>::find_children(
                    node.left_idx,
                    nodes,
                    result,
                    children,
                    octant * 2 );
            }
            if (active1)
            {
                Octree_collector<LEVEL+1,STRIDE>::find_children(
                    node.right_idx,
                    nodes,
                    result,
                    children,
                    octant * 2 + 1 );
            }
        }
    }
};
// Terminal node of Octree_collector's compile-time recursion.
template <uint32_t STRIDE>
struct Octree_collector<3,STRIDE>
{
    static __device__ void find_children(
        const uint32_t        node_index,
        const bintree_node*      nodes,
        Octree_collection*  result,
        uint32_t*             children,
        const uint32_t        octant)
    {
        // we got to one of the octants
        children[ STRIDE * result->node_count ] = node_index;
        result->bitmask |= 1u << octant;
        result->node_count += 1;

        if (nodes[ node_index ].object_idx != 0xFFFFFFFFu) //is_leaf??
            result->leaf_count += 1;
    }
};


// collect octants from a kd-tree
template <uint32_t BLOCK_SIZE>
__global__ void collect_octants_kernel(
    const uint32_t        grid_size,
    const bintree_node*      kd_nodes,
    const uint32_t        in_tasks_count,
    const Split_task*   in_tasks,
    uint32_t*             out_tasks_count,
    Split_task*         out_tasks,
    uint32_t*             out_nodes_count,
    octree_node*        out_nodes,
    const aabb_t*       in_aabb,
    aabb_t*             out_aabb)
{
    const uint32_t LOG_WARP_SIZE = 5;
    const uint32_t WARP_SIZE = 1u << LOG_WARP_SIZE;

    volatile __shared__ uint32_t warp_offset[ BLOCK_SIZE >> LOG_WARP_SIZE ];

    const uint32_t warp_tid = threadIdx.x & (WARP_SIZE-1);
    const uint32_t warp_id  = threadIdx.x >> LOG_WARP_SIZE;

    volatile __shared__ uint32_t sm_red[ BLOCK_SIZE * 2 ];
    volatile uint32_t* warp_red = sm_red + WARP_SIZE * 2 * warp_id;

    __shared__ uint32_t sm_children[ BLOCK_SIZE * 8 ];
    uint32_t* children = sm_children + threadIdx.x;

    // loop through all logical blocks associated to this physical one
    for (uint32_t base_idx = blockIdx.x * BLOCK_SIZE;
        base_idx < in_tasks_count;
        base_idx += grid_size)
    {
        const uint32_t task_id = threadIdx.x + base_idx;

        // uint32_t node;
        Split_task mytask;

        Octree_collection result;
        result.node_count = 0;
        result.leaf_count = 0;
        result.bitmask    = 0;

        // check if the task id is in range, and if so try to collect its treelet
        if (task_id < in_tasks_count)
        {
            // const Split_task in_task = in_tasks[ task_id ];

            mytask = in_tasks[ task_id ];

            Octree_collector<0,BLOCK_SIZE>::find_children(
                mytask.m_input,
                kd_nodes,
                &result,
                children );
        }

        // allocate output nodes, output tasks, and write all leaves
        {
            uint32_t task_count = result.node_count - result.leaf_count;

            uint32_t node_offset = myalloc( result.node_count, out_nodes_count, warp_tid, warp_red, warp_offset + warp_id );
            uint32_t task_offset = myalloc( task_count,        out_tasks_count, warp_tid, warp_red, warp_offset + warp_id );

            // printf("out_task_count: %u\n", *out_tasks_count);

            // write the parent node
            if (task_id < in_tasks_count){
                // out_nodes[ mytask.m_node ] = octree_node( result.bitmask, node_offset, mytask.m_numPts, mytask.m_minIdx, mytask.m_input );
                out_nodes[ mytask.m_node ] = octree_node( result.bitmask, node_offset, mytask.m_numPts, mytask.m_minIdx );
                out_aabb [ mytask.m_node ] = in_aabb[ mytask.m_input ];
            }

            // write out all outputs
            for (uint32_t i = 0; i < result.node_count; ++i)
            {
                const uint32_t  kd_node_index = children[ i * BLOCK_SIZE ];
                const bintree_node kd_node       = kd_nodes[ kd_node_index ];

                if (kd_node.object_idx == 0xFFFFFFFFu) // isInternal??
                    out_tasks[ task_offset++ ] = Split_task( node_offset, kd_node_index, kd_node.numPts, kd_node.min_idx );
                else{ // isLeaf
                    // out_nodes[ node_offset ] = octree_node( kd_node.numPts, kd_node.min_idx, kd_node.object_idx, mytask.m_input );
                    out_nodes[ node_offset ] = octree_node( kd_node.numPts, kd_node.min_idx, kd_node.object_idx );
                    out_aabb [ node_offset ] = in_aabb[ kd_node_index ];
                }

                node_offset++;
            }
        }
    }

// #ifdef DEBUG
//     __syncthreads();
//     if(threadIdx.x == 0)
//         printf("kernel: %u, %u, %u\n out: %u\n", LOG_WARP_SIZE, WARP_SIZE, BLOCK_SIZE, *out_tasks_count);
// #endif
}

// collect octants from a kd-tree
inline void collect_octants(
    const bintree_node* node_list,
    const uint32_t        in_tasks_count,
    const Split_task*   in_tasks,
    uint32_t*             out_tasks_count,
    Split_task*         out_tasks,
    uint32_t*             out_nodes_count,
    octree_node*        out_nodes,
    const aabb_t*       in_aabb,
    aabb_t*             out_aabb)
{
    const uint32_t BLOCK_SIZE = 128;
    int aux;
    cudaOccupancyMaxActiveBlocksPerMultiprocessor( &aux, 
                                                 collect_octants_kernel<BLOCK_SIZE>, BLOCK_SIZE, 
                                                 0);
    const size_t max_blocks = aux;
    const size_t n_blocks   = std::min( max_blocks, size_t((in_tasks_count + BLOCK_SIZE-1) / BLOCK_SIZE) );
    const size_t grid_size  = n_blocks * BLOCK_SIZE;

// #ifdef DEBUG
//         printf("kernel: %zu, %zu, %u\n", grid_size, n_blocks, BLOCK_SIZE);
// #endif

    collect_octants_kernel<BLOCK_SIZE> <<<n_blocks,BLOCK_SIZE>>> (
        grid_size,
        node_list,
        in_tasks_count,
        in_tasks,
        out_tasks_count,
        out_tasks,
        out_nodes_count,
        out_nodes,
        in_aabb,
        out_aabb );

    // cudaThreadSynchronize();
    // cudaDeviceSynchronize();
}


} // namespace octree_builder


Octree_builder::Octree_builder( std::string inputTXT, const uint32_t chunk, cl::sycl::queue& q,
                    const uint32_t no, const uint32_t nn, const uint32_t nin, 
                    const uint32_t nt, const uint32_t non) : 
                    device_queue(q), bintree(inputTXT, chunk, q, no, nn, nin), 
                    b_task_queue1(nt), b_task_queue2(nt),
                    b_counters(3), b_octree(non), b_aabb(non), maxNumTasks(nt), maxNumNodes(non)
{

    // need_space( m_task_queues[0], nt );
    // need_space( m_task_queues[1], nt );

    // m_task_queues[0] = (Split_task*)sycl::malloc_device(nt * sizeof(Split_task), device_queue);
    // m_task_queues[1] = (Split_task*)sycl::malloc_device(nt * sizeof(Split_task), device_queue);

    // need_space( m_counters, 3);
    // need_space( m_octree, numOctreeNodes);
    // need_space( m_aabb, numOctreeNodes );

    reset();

}


void Octree_builder::build()
{


// #ifdef DEBUG
//     std::chrono::time_point<tempo_t> i_start = tempo_t::now();
// #endif

    bintree.build();

    // cudaError_t lastError = cudaGetLastError();
    // if(lastError != cudaSuccess) printf("Error BUILD: %s\n", cudaGetErrorString(lastError));
    // cudaError_t syncError = cudaDeviceSynchronize();
    // if(syncError != cudaSuccess) printf("Error BUILD: %s\n", cudaGetErrorString(syncError));

// #ifdef DEBUG
//     cudaError_t syncError = cudaDeviceSynchronize();
//     if(syncError != cudaSuccess) printf("Error BUILD: %s\n", cudaGetErrorString(syncError));
//     double mytime = cast_t(tempo_t::now() - i_start).count();
//     std::cout << " BINTREE time elapsed: " << mytime << " ms\n";
// #endif

// #ifdef DEBUG
//     i_start = tempo_t::now();
// #endif

    /*Lo primero es reservar espacio para las colas de tareas de entrada y salida,
    y para los contadores*/

    // Split_task* task_queues[2] = {
    //     thrust::raw_pointer_cast( &(m_task_queues[0]).front() ),
    //     thrust::raw_pointer_cast( &(m_task_queues[1]).front() )
    // };

    // sycl::buffer<Split_task, 1> my_tasks[2] = {
    //     sycl::buffer<Split_task, 1>(32),
    //     sycl::buffer<Split_task, 1>(32)
    // };

    sycl::buffer<Split_task, 1> task_queues[2] = {
        b_task_queue1,
        b_task_queue2
    };

    // Split_task* task_queues[2] = {
    //                                 m_task_queues[0],
    //                                 m_task_queues[1]
    // };

    uint32_t in_queue  = 0;
    uint32_t out_queue = 1;

    // uint32_t m_counters_h[3];

    // convert the kd-tree into an octree.
    // m_counters[ in_queue ]  = 1;
    // m_counters[ out_queue ] = 0;
    // m_counters[ 2 ]         = 1; // output node counter

    // bintree_node* node_list_d = (bintree_node*)sycl::malloc_device(bintree.numNodes * sizeof(bintree_node), device_queue);
    // aabb_t* aabb_list_d = (aabb_t*)sycl::malloc_device(bintree.numNodes * sizeof(aabb_t), device_queue);
    
    // bintree_node* node_list_h = (bintree_node*)std::malloc(bintree.numNodes * sizeof(bintree_node));
    // aabb_t* aabb_list_h = (aabb_t*)std::malloc(bintree.numNodes * sizeof(aabb_t));

//     {
//         // auto hnodes = bintree.node_list.get_access<sycl::access::mode::write>();
//         // auto haabbs = bintree.aabb_list.get_access<sycl::access::mode::write>();

//         sycl::host_accessor h_nodes(bintree.node_list, sycl::range<1>(1));
//         sycl::host_accessor h_in_queue(task_queues[ in_queue ], sycl::range<1>(1));

//         h_in_queue[0] = Split_task( 0u, 0u, h_nodes[0].numPts, h_nodes[0].min_idx );

//         // auto task0 = Split_task( 0u, 0u, h_nodes[0].numPts, h_nodes[0].min_idx );

//         // device_queue.memcpy(m_task_queues[ in_queue ], &task0, sizeof(Split_task)).wait();

//         // for(int i=0; i<bintree.numNodes; i++){
//         //     node_list_h[i] = hnodes[i];
//         //     aabb_list_h[i] = haabbs[i];
//         // }

// #ifdef DEBUG
//         printf("Initial Split_task: 0, 0, %u, %u\n", h_nodes[0].numPts, h_nodes[0].min_idx);
// #endif
//     }

    // device_queue.memcpy(node_list_d, node_list_h, bintree.numNodes * sizeof(bintree_node)).wait();
    // device_queue.memcpy(aabb_list_d, aabb_list_h, bintree.numNodes * sizeof(aabb_t)).wait();

    // Split_task init = Split_task( 0u, 0u, node_list_h[0].numPts, node_list_h[0].min_idx );

    // device_queue.memcpy(m_task_queues[ in_queue ], &init, sizeof(Split_task)).wait();

    // m_task_queues[ in_queue ][0] = Split_task( 0u, 0u, node_list_h[0].numPts, node_list_h[0].min_idx );

    // uint32_t level = 0;
    // m_levels[ level++ ] = 0;

    uint32_t counters[3] = {1,0,1};

#ifdef DEBUG
        // printf("Nivel: %u\n", level);
        printf("Tareas in: %u\n", uint32_t(counters[in_queue]));
        printf("Tareas out: %u\n", uint32_t(counters[out_queue]));
        printf("Nodos de salida: %u\n\n", uint32_t(counters[2]));
        // m_levels[ level++ ] = m_counters[2];
#endif

    // loop until there's tasks left in the input queue
    while (counters[ in_queue ])
    {
        // need_space( m_octree, m_counters[2] + m_counters[ in_queue ]*8 );

        // clear the output queue
        // m_counters[ out_queue ] = 0;
        // device_queue.memcpy(m_counters, m_counters_h, 3*sizeof(uint32_t)).wait();


        // octree_builder::collect_octants(
        //     node_list_d,
        //     m_counters_h[ in_queue ],
        //     m_task_queues[ in_queue ],
        //     m_counters + out_queue,
        //     m_task_queues[ out_queue ],
        //     m_counters + 2,
        //     m_octree,
        //     aabb_list_d,
        //     m_aabb );


        const uint32_t BLOCK_SIZE = 128;
        int aux;
        cudaOccupancyMaxActiveBlocksPerMultiprocessor( &aux, 
                                                    octree_builder::collect_octants_kernel<BLOCK_SIZE>, BLOCK_SIZE, 
                                                    0);
        const size_t max_blocks = aux;
        const size_t n_blocks   = std::min( max_blocks, size_t((counters[ in_queue ] + BLOCK_SIZE-1) / BLOCK_SIZE) );
        const size_t grid_size  = n_blocks * BLOCK_SIZE;

// #ifdef DEBUG
//         printf("kernel: %zu, %zu, %u\n", grid_size, n_blocks, BLOCK_SIZE);
// #endif

        // octree_builder::collect_octants_kernel<BLOCK_SIZE> <<<n_blocks,BLOCK_SIZE>>> (
        //     grid_size,
        //     node_list_d,
        //     m_counters_h[ in_queue ],
        //     m_task_queues[ in_queue ],
        //     m_counters + out_queue,
        //     m_task_queues[ out_queue ],
        //     m_counters + 2,
        //     m_octree,
        //     aabb_list_d,
        //     m_aabb );


        device_queue.submit([&](sycl::handler& h) {

            // auto accN = bintree.node_list.get_access<sycl::access::mode::read>(h);
            // auto accA = bintree.aabb_list.get_access<sycl::access::mode::read>(h);
            sycl::accessor accN(bintree.node_list, h, sycl::read_only);
            sycl::accessor accAB(bintree.aabb_list, h, sycl::read_only);
            sycl::accessor accC(b_counters, h, sycl::write_only);
            sycl::accessor accO(b_octree, h, sycl::write_only);
            sycl::accessor accAO(b_aabb, h, sycl::write_only);
            sycl::accessor accIN(task_queues[ in_queue ], h, sycl::read_only);
            sycl::accessor accOUT(task_queues[ out_queue ], h, sycl::write_only);

            h.interop_task([=](sycl::interop_handler ih) {

                auto dN = reinterpret_cast<node_t*>(ih.get_mem<sycl::backend::cuda>(accN));
                auto dAB = reinterpret_cast<aabb_t*>(ih.get_mem<sycl::backend::cuda>(accAB));
                auto dC = reinterpret_cast<uint32_t*>(ih.get_mem<sycl::backend::cuda>(accC));
                auto dO = reinterpret_cast<octree_node*>(ih.get_mem<sycl::backend::cuda>(accO));
                auto dAO = reinterpret_cast<aabb_t*>(ih.get_mem<sycl::backend::cuda>(accAO));
                auto dIN = reinterpret_cast<Split_task*>(ih.get_mem<sycl::backend::cuda>(accIN));
                auto dOUT = reinterpret_cast<Split_task*>(ih.get_mem<sycl::backend::cuda>(accOUT));

                octree_builder::collect_octants_kernel<BLOCK_SIZE> <<<n_blocks,BLOCK_SIZE>>> (
                    grid_size,
                    dN,
                    counters[ in_queue ],
                    dIN,
                    dC + out_queue,
                    dOUT,
                    dC + 2,
                    dO,
                    dAB,
                    dAO );

                // octree_builder::collect_octants(
                //     dN,
                //     m_counters_h[ in_queue ],
                //     m_task_queues[ in_queue ],
                //     m_counters + out_queue,
                //     m_task_queues[ out_queue ],
                //     m_counters + 2,
                //     m_octree,
                //     dA,
                //     m_aabb );

                // octree_builder::collect_octants(
                //     dN,
                //     m_counters[ in_queue ],
                //     task_queues[ in_queue ],
                //     thrust::raw_pointer_cast( &m_counters.front() ) + out_queue,
                //     task_queues[ out_queue ],
                //     thrust::raw_pointer_cast( &m_counters.front() ) + 2,
                //     thrust::raw_pointer_cast( &m_octree.front() ),
                //     dA,
                //     thrust::raw_pointer_cast( &m_aabb.front() ) );

            });

        });

        // sycl::host_accessor hnodes(bintree.node_list);
        // sycl::host_accessor haabbs(bintree.aabb_list);

        // cudaDeviceSynchronize();

        // cudaError_t syncError = cudaDeviceSynchronize();
        // if(syncError != cudaSuccess) printf("Error BUILD: %s\n", cudaGetErrorString(syncError));

#ifdef DEBUG
        cudaError_t lastError = cudaGetLastError();
        if(lastError != cudaSuccess) printf("Error BUILD: %s\n", cudaGetErrorString(lastError));
        cudaError_t syncError = cudaDeviceSynchronize();
        if(syncError != cudaSuccess) printf("Error BUILD: %s\n", cudaGetErrorString(syncError));

        // printf("Nivel: %u\n", level);
        // printf("Tareas in: %u\n", uint32_t(m_counters[in_queue]));
        // printf("Tareas out: %u\n", uint32_t(m_counters[out_queue]));
        // printf("Nodos de salida; %u\n\n", uint32_t(m_counters[2]));
        // m_levels[ level++ ] = m_counters[2];
#endif

        // device_queue.memcpy(m_counters_h, m_counters, 3*sizeof(uint32_t)).wait();

        // printf("En el nivel %u he encontrado %u leafs\n", level, m_counters[2] - m_counters[out_queue]);

        // swap the input and output queues
        std::swap( in_queue, out_queue );

        {
            cudaDeviceSynchronize();
            
            sycl::host_accessor h_counters(b_counters, sycl::read_write);

            counters[ in_queue ] = h_counters[ in_queue ]; // in_queue            
            // counters[ out_queue ] = h_counters[ out_queue ]; // out_queue            
            counters[ 2 ] = h_counters[ 2 ]; 

            h_counters[ out_queue ] = 0u; // out_queue
        }

#ifdef DEBUG
        // printf("Nivel: %u\n", level);
        printf("Tareas in: %u\n", uint32_t(counters[in_queue]));
        // printf("Tareas out: %u\n", uint32_t(counters[out_queue]));
        printf("Nodos de salida: %u\n\n", uint32_t(counters[2]));
        // m_levels[ level++ ] = m_counters[2];
#endif

    }

    // cudaError_t syncError = cudaDeviceSynchronize();
    // if(syncError != cudaSuccess) printf("Error BUILD: %s\n", cudaGetErrorString(syncError));

    // {
    //     /* Inicializo aquí raíz porque me da lo mismo, solo me importan los dos primeros
    //     campos de Split_task que SEGURO son 0 y 0*/
        
    //     sycl::host_accessor h_nodes(bintree.node_list, sycl::range<1>(1));
    //     sycl::host_accessor h_octree(b_octree, sycl::range<1>(1));

    //     h_octree[0].numPts = h_nodes[0].numPts;
    //     h_octree[0].min_idx = h_nodes[0].min_idx;

    // }


    m_leaf_count = bintree.numLeafs;
    m_node_count = counters[2];

    // for (; level < 64; ++level)
    //     m_levels[ level ] = m_node_count;

// #ifdef DEBUG
//     double mytime = cast_t(tempo_t::now() - i_start).count();
//     std::cout << " OCTREE time elapsed: " << mytime << " ms\n";
// #endif

    // free(node_list_d, device_queue);
    // free(aabb_list_d, device_queue);
    // free(node_list_h);
    // free(aabb_list_h);

    return;
}

void Octree_builder::reset()
{

    bintree.reset();
    
    {
        sycl::host_accessor h_counters(b_counters);
        h_counters[ 0 ]  = 1; // in_queue
        h_counters[ 1 ]  = 0; // out_queue
        h_counters[ 2 ]  = 1; // output node counter
    }

    device_queue.submit([&](sycl::handler& h) {

        sycl::accessor accIN(b_task_queue1, h, sycl::write_only);
        sycl::accessor accOUT(b_task_queue2, h, sycl::write_only);

        h.interop_task([=](sycl::interop_handler ih) {

            auto dIN = reinterpret_cast<Split_task*>(ih.get_mem<sycl::backend::cuda>(accIN));
            auto dOUT = reinterpret_cast<Split_task*>(ih.get_mem<sycl::backend::cuda>(accOUT));

            init_tasks<<<30,128>>>(dIN, dOUT, maxNumTasks);

        });
    });

    device_queue.submit([&](sycl::handler& h) {

        sycl::accessor accN(b_octree, h, sycl::write_only);
        sycl::accessor accA(b_aabb, h, sycl::write_only);

        h.interop_task([=](sycl::interop_handler ih) {

            auto dN = reinterpret_cast<octree_node*>(ih.get_mem<sycl::backend::cuda>(accN));
            auto dA = reinterpret_cast<aabb_t*>(ih.get_mem<sycl::backend::cuda>(accA));

            init_octree<<<30,128>>>(dN, dA, maxNumNodes, default_aabb);

        });
    });

    // device_queue.memset(m_task_queues[0], 0u, maxNumTasks * sizeof(Split_task)).wait();
    // device_queue.memset(m_task_queues[1], 0u, maxNumTasks * sizeof(Split_task)).wait();

    cudaError_t lastError = cudaGetLastError();
    if(lastError != cudaSuccess) printf("Error RESET: %s\n", cudaGetErrorString(lastError));
    cudaError_t syncError = cudaDeviceSynchronize();
    if(syncError != cudaSuccess) printf("Error RESET: %s\n", cudaGetErrorString(syncError));
}
