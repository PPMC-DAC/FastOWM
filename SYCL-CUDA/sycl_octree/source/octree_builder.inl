#include <algorithm>

/// intra-warp inclusive scan
///
/// \param val      per-threrad input value
/// \param tidx     warp thread index
/// \param red      scan result storage (2*WARP_SIZE elements)
template <typename T> 
inline T scan_warp(T val, const int32_t tidx, volatile T *red)
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

template <typename T> 
inline T scan_warp(T val, const int32_t tidx, volatile T *red, sycl::nd_item<1>& it )
{
    // pad initial segment with zeros
    red[tidx] = 0;
#ifdef CPU
    red += 8;
#else
    red += 32;
#endif
    // Hillis-Steele scan
    red[tidx] = val;
    it.barrier(sycl::access::fence_space::local_space);

#ifdef CPU
    for(int offset = 1; offset<8; offset *= 2)
#else
    for(int offset = 1; offset<32; offset *= 2)
#endif
    {
        val += red[tidx-offset];  red[tidx] = val;
        it.barrier(sycl::access::fence_space::local_space);
    }

	return val;
}


/// return the total from a scan_warp
///
/// \param red      scan result storage
template <typename T> 
inline T scan_warp_total(volatile T *red) { 
#ifdef CPU
    return red[15]; 
#else
    return red[63];
#endif
    }


#ifdef NVIDIA
/// alloc n elements per thread from a common pool, using a synchronous warp scan
///
/// \param n                number of elements to alloc
/// \param warp_tid         warp thread index
/// \param warp_red         temporary warp scan storage (2*WARP_SIZE elements)
/// \param warp_broadcast   temporary warp broadcasting storage
inline
uint32_t myalloc(uint32_t n, uint32_t* pool, const int32_t warp_tid, volatile uint32_t* warp_red, volatile uint32_t* warp_broadcast)
{
    uint32_t warp_scan  = scan_warp( n, warp_tid, warp_red ) - n;
    uint32_t warp_count = scan_warp_total( warp_red );
    if (warp_tid == 0)
        *warp_broadcast = atomic_acc(*pool).fetch_add(warp_count);

    return *warp_broadcast + warp_scan;
}
#else
inline
uint32_t myalloc(uint32_t n, uint32_t* pool, const int32_t warp_tid, volatile uint32_t* warp_red, volatile uint32_t* warp_broadcast, sycl::nd_item<1>& it)
{
    uint32_t warp_scan  = scan_warp( n, warp_tid, warp_red, it ) - n;

    uint32_t warp_count = scan_warp_total( warp_red );

    if (warp_tid == 0)
        // *warp_broadcast = atomic_wg(*pool).fetch_add(warp_count);
        *warp_broadcast = sycl::atomic<uint32_t>(sycl::global_ptr<uint32_t>(pool)).fetch_add(warp_count);

    it.barrier(sycl::access::fence_space::local_space);

    return *warp_broadcast + warp_scan;
}
#endif

// template <typename T> 
// inline T _scan(T val, const int32_t tidx, volatile T *red, sycl::nd_item<1>& it )
// {
//     int pout = 0, pin = 1;
//     red[tidx] = val;
//     it.barrier(sycl::access::fence_space::local_space);

//     red[pout*32 + tidx] = (tidx>0)? red[tidx-1] : 0;
//     it.barrier(sycl::access::fence_space::local_space);

//     for(int offset=1; offset<32; offset*=2)
//     {
//         pout = 1-pout;
//         pin = 1-pout;

//         if(tidx >= offset)
//             red[pout*32+tidx] += red[pin*32+tidx - offset];
//         else
//             red[pout*32+tidx] = red[pin*32+tidx];

//         it.barrier(sycl::access::fence_space::local_space);
//     }
    
//     return red[pout*32+tidx];

// }





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
    static void find_children(
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
    static void find_children(
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
inline void collect_octants(
    const bintree_node* node_list,
    const uint32_t        in_tasks_count,
    const Split_task*   in_tasks,
    uint32_t*             out_tasks_count,
    Split_task*         out_tasks,
    uint32_t*             out_nodes_count,
    octree_node*        out_nodes,
    const aabb_t*       in_aabb,
    aabb_t*             out_aabb,
    sycl::queue&        device_queue)
{
#ifdef CPU
    const uint32_t LOG_WARP_SIZE = 3;
#else
    const uint32_t LOG_WARP_SIZE = 5;
#endif
    const uint32_t WARP_SIZE = 1u << LOG_WARP_SIZE;

    const uint32_t BLOCK_SIZE = 128;
    // const size_t n_blocks   = 1;
    const size_t grid_size  = BLOCK_SIZE;

    device_queue.submit([&](sycl::handler& cgh) {

        sycl::local_accessor<uint32_t, 1>
            warp_offset(sycl::range<1>(BLOCK_SIZE >> LOG_WARP_SIZE), cgh);

        sycl::local_accessor<uint32_t, 1>
            sm_red(sycl::range<1>(BLOCK_SIZE*2), cgh);

        sycl::local_accessor<uint32_t, 1>
            sm_children(sycl::range<1>(BLOCK_SIZE*8), cgh);

        // sycl::stream out(16 * 1024, 16 * 1024, cgh);

        cgh.parallel_for<class collect_kernel>(
            sycl::nd_range<1>(grid_size, BLOCK_SIZE),
            [=] (sycl::nd_item<1> it) {

                const uint32_t warp_tid = it.get_local_id(0) & (WARP_SIZE-1);
                const uint32_t warp_id  = it.get_local_id(0) >> LOG_WARP_SIZE;

                uint32_t* warp_red = &(sm_red[0]) + WARP_SIZE * 2 * warp_id;

                // if( it.get_local_id(0) == 0 ) out << "THREAD 0 speaking: \n";

                uint32_t* children = &(sm_children[0]) + it.get_local_id(0);

                // loop through all logical blocks associated to this physical one
                for (uint32_t base_idx = it.get_group(0) * BLOCK_SIZE;
                    base_idx < in_tasks_count;
                    base_idx += grid_size)
                {
                    const uint32_t task_id = it.get_local_id(0) + base_idx;

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
                            node_list,
                            &result,
                            children );
                    }

                    // allocate output nodes, output tasks, and write all leaves
                    {
                        uint32_t task_count = result.node_count - result.leaf_count;

                        // if( it.get_local_id(0) == 0 ) out << "   " << task_count << " new tasks\n";
                        // out << "Thread " << task_id << ": " << task_count << " new tasks\n";

#ifdef NVIDIA
                        uint32_t node_offset = myalloc( result.node_count, out_nodes_count, warp_tid, warp_red, &(warp_offset[0]) + warp_id );
                        uint32_t task_offset = myalloc( task_count, out_tasks_count, warp_tid, warp_red, &(warp_offset[0]) + warp_id );
#else
                        uint32_t node_offset = myalloc( result.node_count, out_nodes_count, warp_tid, warp_red, &(warp_offset[0]) + warp_id, it );
                        uint32_t task_offset = myalloc( task_count, out_tasks_count, warp_tid, warp_red, &(warp_offset[0]) + warp_id, it );
#endif

                        // if( it.get_local_id(0) == 0 ) out << "   " << *out_nodes_count << " nodes saved\n" 
                                                        // << "   " << node_offset << " node offset\n"
                                                        // << "   " << *out_tasks_count << " tasks saved\n"
                                                        // << "   " << task_offset << " task offset\n";

                        // write the parent node
                        if (task_id < in_tasks_count){
                            // out_nodes[ mytask.m_node ] = octree_node( result.bitmask, node_offset, mytask.m_numPts, mytask.m_minIdx, mytask.m_input );
                            out_nodes[ mytask.m_node ] = octree_node( result.bitmask, node_offset, mytask.m_numPts, mytask.m_minIdx );
                            out_aabb [ mytask.m_node ] = in_aabb[ mytask.m_input ];
                        }

                        // write out all outputs
                        for (uint32_t i = 0; i < result.node_count; ++i)
                        {
                            // if( it.get_local_id(0) == 0 ) out << "   node " << i+1 << "\n";

                            const uint32_t  kd_node_index = children[ i * BLOCK_SIZE ];
                            const bintree_node kd_node       = node_list[ kd_node_index ];

                            if (kd_node.object_idx == 0xFFFFFFFFu){ // isInternal??
                                out_tasks[ task_offset++ ] = Split_task( node_offset, kd_node_index, kd_node.numPts, kd_node.min_idx );
                                // out << "Thread " << task_id << ": " << task_offset << " new tasks\n";
                            }else{ // isLeaf
                                // out_nodes[ node_offset ] = octree_node( kd_node.numPts, kd_node.min_idx, kd_node.object_idx, mytask.m_input );
                                out_nodes[ node_offset ] = octree_node( kd_node.numPts, kd_node.min_idx, kd_node.object_idx );
                                out_aabb [ node_offset ] = in_aabb[ kd_node_index ];
                            }

                            node_offset++;
                        }

                        // if( it.get_local_id(0) == 0 ) out << "   " << task_offset << " task offset\n";
                        // out << "Thread " << task_id << ": " << task_offset << " new tasks\n";
                    }
                    
                } //for
                
                return;
            }
        );
    
    }).wait_and_throw();

}


} // namespace octree_builder


Octree_builder::Octree_builder( std::string inputTXT, const uint32_t chunk, sycl::queue& q) : 
                bintree(inputTXT, chunk, q), device_queue(q)
    {
        // Numero minimo de nodos en bloques de 8 que contienen a todas las hojas
        uint32_t aux = uint32_t((bintree.numLeafs-1)/8) + 1u;
        // pero necesito todos los nodos intermedios también
        int i=0;
        uint32_t num_tasks = 0;
        uint32_t pot = 0;
        for(; pot<aux; i++) {
            pot = uint32_t(pow(8,i));
            num_tasks += pot;
        }
        // Procesar un nodo hoja no es una tarea para nosotros
        // need_space( m_task_queues[0], num_tasks );
        // need_space( m_task_queues[1], num_tasks );
        maxNumTasks = num_tasks;
        m_task_queues[0] = (Split_task*)sycl_builder::mallocWrap(maxNumTasks*sizeof(Split_task), device_queue);
        m_task_queues[1] = (Split_task*)sycl_builder::mallocWrap(maxNumTasks*sizeof(Split_task), device_queue);
        
        // Pero un nodo hoja si va a ser un nuevo nodo del octree
        // num_tasks += bintree.numLeafs;
        maxNumNodes = maxNumTasks + bintree.numLeafs;

    #ifdef DEBUG
        printf("Necesito %d niveles en el arbol\n", i);
        printf("Voy a tener como máximo %u nodos\n", num_tasks);
    #endif

        // need_space( m_counters, 3);
        // need_space( m_octree, num_tasks);
        // need_space( m_aabb, num_tasks );
        m_counters = (uint32_t*)sycl_builder::mallocWrap(3*sizeof(uint32_t), device_queue);
        m_octree = (octree_node*)sycl_builder::mallocWrap(maxNumNodes*sizeof(octree_node), device_queue);
        m_aabb = (aabb_t*)sycl_builder::mallocWrap(maxNumNodes*sizeof(aabb_t), device_queue);

        reset();
    }



void Octree_builder::build()
{

    bintree.build();

#ifdef DEBUG
    std::chrono::time_point<tempo_t> i_start = tempo_t::now();
#endif

    Split_task* task_queues[2] = {
        m_task_queues[0],
        m_task_queues[1]
    };

    uint32_t in_queue  = 0;
    uint32_t out_queue = 1;

    // convert the kd-tree into an octree.
#ifdef DEVICE
    uint32_t counters[3] = {1,0,1};
#endif

    uint32_t num_internal_nodes = 1;

    // m_task_queues[ in_queue ][0] = Split_task( 0u, 0u, bintree.node_list[0].numPts, bintree.node_list[0].min_idx );

    // uint32_t level = 0;
    // m_levels[ level++ ] = 0;

#ifdef DEBUG

#ifdef DEVICE
    printf("Tareas in: %u\n", uint32_t(counters[in_queue]));
    printf("Tareas out: %u\n", uint32_t(counters[out_queue]));
    printf("Nodos de salida: %u\n\n", uint32_t(counters[2]));
#else
    printf("Tareas in: %u\n", uint32_t(m_counters[in_queue]));
    printf("Tareas out: %u\n", uint32_t(m_counters[out_queue]));
    printf("Nodos de salida: %u\n\n", uint32_t(m_counters[2]));
#endif

#endif

    // loop until there's tasks left in the input queue
#ifdef DEVICE
    while (counters[ in_queue ])
#else
    while (m_counters[ in_queue ])
#endif
    {

        // need_space( m_octree, m_counters[2] + m_counters[ in_queue ]*8 );

        // clear the output queue
        // counters[ out_queue ] = 0;

        octree_builder::collect_octants(
            bintree.node_list,
#ifdef DEVICE
            counters[ in_queue ],
#else
            m_counters[ in_queue ],
#endif
            task_queues[ in_queue ],
            m_counters + out_queue,
            task_queues[ out_queue ],
            m_counters + 2,
            m_octree,
            bintree.aabb_list,
            m_aabb,
            device_queue
        );

        // printf("En el nivel %u he encontrado %u leafs\n", level, m_counters[2] - m_counters[out_queue]);

        // swap the input and output queues
        std::swap( in_queue, out_queue );

#ifdef DEVICE
        device_queue.memcpy(counters, m_counters, 3*sizeof(uint32_t)).wait();
#endif

#ifdef DEBUG

#ifdef DEVICE
    printf("Tareas in: %u\n", uint32_t(counters[in_queue]));
    printf("Tareas out: %u\n", uint32_t(counters[out_queue]));
    printf("Nodos de salida: %u\n\n", uint32_t(counters[2]));
#else
    printf("Tareas in: %u\n", uint32_t(m_counters[in_queue]));
    printf("Tareas out: %u\n", uint32_t(m_counters[out_queue]));
    printf("Nodos de salida: %u\n\n", uint32_t(m_counters[2]));
#endif

#endif

#ifdef DEVICE
        num_internal_nodes += counters[in_queue];

        counters[ out_queue ] = 0;

        device_queue.memcpy(m_counters, counters, 3*sizeof(uint32_t)).wait();
#else
        num_internal_nodes += m_counters[ in_queue ];

        m_counters[ out_queue ] = 0;
#endif

    }

#ifdef DEVICE
    m_leaf_count = counters[2] - num_internal_nodes;
    m_node_count = counters[2];
#else
    m_leaf_count = m_counters[2] - num_internal_nodes;
    m_node_count = m_counters[2];
#endif

#ifdef DEBUG
    double mytime = cast_t(tempo_t::now() - i_start).count();
    std::cout << " OCTREE time elapsed: " << mytime << " ms\n";
#endif

    return;
}

class _reset_tasks{
  
  public:
    _reset_tasks( Split_task* q1, Split_task* q2) : 
        in_queue(q1), out_queue(q2) {}

    void operator()(sycl::id<1> idx) const
    {
      in_queue[idx] = Split_task(0u,0u,0u,0u);
      out_queue[idx] = Split_task(0u,0u,0u,0u);

      return;
    }
  
  private:
    Split_task*  in_queue;
    Split_task*  out_queue;

};

class _reset_octree{
  
  public:
    _reset_octree( octree_node* n, aabb_t* bb, const aabb_t def_bb) : 
        octree_list(n), aabb_list(bb), default_aabb(def_bb) {}

    void operator()(sycl::id<1> idx) const
    {
      octree_list[idx] = octree_node(0u,0u,0u,0u);
      aabb_list[idx] = default_aabb;

      return;
    }
  
  private:
    octree_node*  octree_list;
    aabb_t*       aabb_list;
    const aabb_t default_aabb;

};

void Octree_builder::reset()
{
    bintree.reset();

#ifdef DEVICE
    uint32_t counters[3] = {1,0,1};

    device_queue.memcpy(m_counters, counters, 3*sizeof(uint32_t)).wait();
#else
    m_counters[ 0 ]     = 1;
    m_counters[ 1 ]     = 0;
    m_counters[ 2 ]     = 1; // output node counter
#endif

    device_queue.submit([&](sycl::handler &cgh) {

        cgh.parallel_for(sycl::range<1>(maxNumTasks), 
                            _reset_tasks(m_task_queues[0], m_task_queues[1]));

    });

    device_queue.submit([&](sycl::handler &cgh) {

        cgh.parallel_for(sycl::range<1>(maxNumNodes), 
                            _reset_octree(m_octree, m_aabb, default_aabb));

    });

    device_queue.wait_and_throw();

}
