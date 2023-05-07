#include <algorithm>
#include <cub/cub.cuh>
#include <cub/device/device_radix_sort.cuh>

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
__device__ __forceinline__ constexpr T scan_warp_total(volatile T *red) { return red[63]; }

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
    static __host__ __device__ void find_children(
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
            // result->bitmask |= 1u << (octant << (3 - LEVEL));
            /*Puedo comprobar comprobar cual es el último activado y ponerlo en el siguiente*/
            int j=0;
            while(result->bitmask & (1u << j)) j++;
            result->bitmask |= 1u << j;
            
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
    static __host__ __device__ void find_children(
        const uint32_t        node_index,
        const bintree_node*      nodes,
        Octree_collection*  result,
        uint32_t*             children,
        const uint32_t        octant)
    {
        // we got to one of the octants
        children[ STRIDE * result->node_count ] = node_index;

        /*En vez de hacer esto necesito una forma de que los octantes estén correlativos
        y así poder acceder posteriormente más fácil*/
        // result->bitmask |= 1u << octant;

        /*Puedo comprobar comprobar cual es el último activado y ponerlo en el siguiente*/
        int j=0;
        while(result->bitmask & (1u << j)) j++;
        result->bitmask |= 1u << j;

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
    cudaDeviceSynchronize();
}


} // namespace octree_builder


Octree_builder::Octree_builder( std::string inputTXT, const uint32_t chunk) : bintree(inputTXT, chunk)
{
    // Numero minimo de nodos en bloques de 8 que contienen a todas las hojas
    uint32_t aux = uint32_t((bintree.numLeafs-1)/8) + 1u;
    // pero necesito todos los nodos intermedios también
    int i=0;
    maxNumTasks = 0;
    uint32_t pot = 0;
    for(; pot<aux; i++) {
        pot = uint32_t(pow(8,i));
        maxNumTasks += pot;
    }
    // Procesar un nodo hoja no es una tarea para nosotros
    need_space( m_task_queues[0], maxNumTasks );
    need_space( m_task_queues[1], maxNumTasks );
    // cudaMalloc(&(b_task_queues[0]), maxNumTasks * sizeof(Split_task));
    // cudaMalloc(&(b_task_queues[1]), maxNumTasks * sizeof(Split_task));

// #ifdef DEBUG
//     printf("Necesito %d niveles en el arbol\n", i);
//     printf("Voy a tener como máximo %u nodos\n", num_nodes);
// #endif

    need_space( m_counters, 3);
    // need_space( m_octree, num_nodes);
    // need_space( m_aabb, num_nodes );
    // cudaMallocManaged(&b_counters, 3*sizeof(uint32_t));

    // Pero un nodo hoja sí que va a ser un nuevo nodo del octree
    maxNumNodes = maxNumTasks + bintree.numLeafs;
    cudaMallocManaged(&b_octree, maxNumNodes*sizeof(octree_node));
    cudaMallocManaged(&b_aabb, maxNumNodes*sizeof(aabb_t));

    reset();

}


void Octree_builder::build()
{

    bintree.build();

#ifdef DEBUG
    std::chrono::time_point<tempo_t> i_start = tempo_t::now();
#endif

    /*Lo primero es reservar espacio para las colas de tareas de entrada y salida,
    y para los contadores*/
    // thrust::device_vector<Split_task> m_task_queues[2];
    // thrust::device_vector<uint32_t>     m_counters(3);
    // thrust::device_vector<octree_node>      m_octree;

    Split_task* task_queues[2] = {
        thrust::raw_pointer_cast( &(m_task_queues[0]).front() ),
        thrust::raw_pointer_cast( &(m_task_queues[1]).front() )
    };
    // Split_task* task_queues[2] = {
    //     b_task_queues[0],
    //     b_task_queues[1]
    // };

    uint32_t in_queue  = 0;
    uint32_t out_queue = 1;

    // convert the kd-tree into an octree.
    // m_counters[ in_queue ]  = 1;
    // m_counters[ out_queue ] = 0;
    // m_counters[ 2 ]         = 1; // output node counter

    uint32_t num_internal_nodes = 1;

    m_task_queues[ in_queue ][0] = Split_task( 0u, 0u, bintree.node_list[0].numPts, bintree.node_list[0].min_idx );

    // uint32_t level = 0;
    // m_levels[ level++ ] = 0;

    // loop until there's tasks left in the input queue
    while (m_counters[ in_queue ])
    {

#ifdef DEBUG
        // printf("Nivel: %u\n", level);
        printf("Tareas IN: %u\n", uint32_t(m_counters[in_queue]));
        printf("Tareas OUT: %u\n", uint32_t(m_counters[out_queue]));
        printf("Nodos de salida; %u\n", uint32_t(m_counters[2]));
        // m_levels[ level++ ] = m_counters[2];
#endif

        // need_space( m_octree, m_counters[2] + m_counters[ in_queue ]*8 );

        octree_builder::collect_octants(
            bintree.node_list,
            m_counters[ in_queue ],
            task_queues[ in_queue ],
            thrust::raw_pointer_cast( &(m_counters).front() ) + out_queue,
            task_queues[ out_queue ],
            thrust::raw_pointer_cast( &(m_counters).front() ) + 2,
            b_octree,
            bintree.aabb_list,
            b_aabb
        );

        // cudaMemPrefetchAsync(m_counters, 3*sizeof(uint32_t), cudaCpuDeviceId);

        // printf("En el nivel %u he encontrado %u leafs\n", level, m_counters[2] - m_counters[out_queue]);

        // swap the input and output queues
        std::swap( in_queue, out_queue );

        num_internal_nodes += m_counters[ in_queue ];

        // clear the output queue
        m_counters[ out_queue ] = 0;

        // cudaMemPrefetchAsync(m_counters, 3*sizeof(uint32_t), 0);


    }

#ifdef DEBUG
    // printf("Nivel: %u\n", level);
    printf("Tareas IN: %u\n", uint32_t(m_counters[ in_queue ]));
    printf("Tareas OUT: %u\n", uint32_t(m_counters[ out_queue ]));
    printf("Nodos de salida; %u\n", uint32_t(m_counters[2]));
    // m_levels[ level++ ] = m_counters[2];
#endif

    m_leaf_count = m_counters[2] - num_internal_nodes;
    m_node_count = m_counters[2];

    // for (; level < 64; ++level)
    //     m_levels[ level ] = m_node_count;

#ifdef DEBUG
    double mytime = cast_t(tempo_t::now() - i_start).count();
    std::cout << " OCTREE time elapsed: " << mytime << " ms\n";
#endif

    return;
}

void Octree_builder::reset()
{
  
    bintree.reset();
  
    m_counters[ 0 ]  = 1; // in_queue
    m_counters[ 1 ]  = 0; // out_queue
    m_counters[ 2 ]  = 1; // output node counter

    init_tasks<<<30,128>>>( thrust::raw_pointer_cast( &(m_task_queues[0]).front() ), 
                            thrust::raw_pointer_cast( &(m_task_queues[1]).front() ),
                            maxNumTasks);

    init_octree<<<30,128>>>(b_octree, b_aabb, maxNumNodes, default_aabb);

    cudaError_t lastError = cudaGetLastError();
    if(lastError != cudaSuccess) printf("Error RESET: %s\n", cudaGetErrorString(lastError));
    cudaError_t syncError = cudaDeviceSynchronize();
    if(syncError != cudaSuccess) printf("Error RESET: %s\n", cudaGetErrorString(syncError));

}
