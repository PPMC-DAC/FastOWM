#ifndef LBVH_BVH_CUH
#define LBVH_BVH_CUH
#include "aabb.cuh"
#include "morton_code.cuh"
#include <thrust/swap.h>
#include <thrust/pair.h>
#include <thrust/tuple.h>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/functional.h>
#include <thrust/scan.h>
#include <thrust/sort.h>
#include <thrust/fill.h>
#include <thrust/for_each.h>
#include <thrust/transform.h>
#include <thrust/reduce.h>
#include <thrust/iterator/constant_iterator.h>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/execution_policy.h>
#include <chrono>

using tempo_t = std::chrono::steady_clock;

using cast_t = std::chrono::duration<double, std::milli>;


namespace lbvh
{
namespace detail
{
struct node
{
    std::uint32_t parent_idx; // parent node
    std::uint32_t left_idx;   // index of left  child node
    std::uint32_t right_idx;  // index of right child node
    std::uint32_t object_idx; // == 0xFFFFFFFF if internal node.
};

// a set of pointers to use it on device.
template<typename Real, typename Object, bool IsConst>
struct basic_device_bvh;
template<typename Real, typename Object>
struct basic_device_bvh<Real, Object, false>
{
    using real_type  = Real;
    using aabb_type  = aabb<real_type>;
    using node_type  = detail::node;
    using index_type = std::uint32_t;
    using object_type = Object;

    uint32_t num_nodes;   // (# of internal node) + (# of leaves), 2N+1
    uint32_t num_objects; // (# of leaves), the same as the number of objects

    node_type *  nodes;
    aabb_type *  aabbs;
    object_type* objects;
};
template<typename Real, typename Object>
struct basic_device_bvh<Real, Object, true>
{
    using real_type  = Real;
    using aabb_type  = aabb<real_type>;
    using node_type  = detail::node;
    using index_type = std::uint32_t;
    using object_type = Object;

    uint32_t num_nodes;  // (# of internal node) + (# of leaves), 2N+1
    uint32_t num_objects;// (# of leaves), the same as the number of objects

    node_type   const* nodes;
    aabb_type   const* aabbs;
    object_type const* objects;
};

template<typename UInt>
__device__ __forceinline__
uint2 determine_range(UInt const* node_code,
        const uint32_t num_leaves, uint32_t idx)
{
    if(idx == 0)
    {
        return make_uint2(0, num_leaves-1);
    }

    // determine direction of the range
    const UInt self_code = node_code[idx];
    const int L_delta = common_upper_bits(self_code, node_code[idx-1]);
    const int R_delta = common_upper_bits(self_code, node_code[idx+1]);
    const int d = (R_delta > L_delta) ? 1 : -1;

    // Compute upper bound for the length of the range

    const int delta_min = thrust::min(L_delta, R_delta);
    int l_max = 2;
    int delta = -1;
    int i_tmp = idx + d * l_max;
    if(0 <= i_tmp && i_tmp < num_leaves)
    {
        delta = common_upper_bits(self_code, node_code[i_tmp]);
    }
    while(delta > delta_min)
    {
        l_max <<= 1;
        i_tmp = idx + d * l_max;
        delta = -1;
        if(0 <= i_tmp && i_tmp < num_leaves)
        {
            delta = common_upper_bits(self_code, node_code[i_tmp]);
        }
    }

    // Find the other end by binary search
    int l = 0;
    int t = l_max >> 1;
    while(t > 0)
    {
        i_tmp = idx + (l + t) * d;
        delta = -1;
        if(0 <= i_tmp && i_tmp < num_leaves)
        {
            delta = common_upper_bits(self_code, node_code[i_tmp]);
        }
        if(delta > delta_min)
        {
            l += t;
        }
        t >>= 1;
    }
    uint32_t jdx = idx + l * d;
    if(d < 0)
    {
        thrust::swap(idx, jdx); // make it sure that idx < jdx
    }
    return make_uint2(idx, jdx);
}

template<typename UInt>
__device__ __forceinline__
uint32_t find_split(UInt const* node_code, const uint32_t num_leaves,
    const uint32_t first, const uint32_t last) noexcept
{
    const UInt first_code = node_code[first];
    const UInt last_code  = node_code[last];
    if (first_code == last_code)
    {
        return (first + last) >> 1;
    }   
    const int delta_node = common_upper_bits(first_code, last_code);
// #ifdef DEBUG
//     printf("%u, %u delta: %d\n",first, last, delta_node);
// #endif

    // binary search...
    int split  = first;
    int stride = last - first;
    do
    {
        stride = (stride + 1) >> 1;
        const int middle = split + stride;
        if (middle < last)
        {
            const int delta = common_upper_bits(first_code, node_code[middle]);
            if (delta > delta_node)
            {
                split = middle;
            }
        }
    }
    while(stride > 1);

    return split;
}
template<typename Real, typename Object, bool IsConst, typename UInt>
void construct_internal_nodes(const basic_device_bvh<Real, Object, IsConst>& self,
        UInt const* node_code, const uint32_t num_objects)
{
    thrust::for_each(thrust::device,
        thrust::make_counting_iterator<uint32_t>(0),
        thrust::make_counting_iterator<uint32_t>(num_objects - 1),
        [self, node_code, num_objects] __device__ (const uint32_t idx)
        {
            self.nodes[idx].object_idx = 0xFFFFFFFF; //  internal nodes

            const uint2 ij  = determine_range(node_code, num_objects, idx);
            const int gamma = find_split(node_code, num_objects, ij.x, ij.y);

            self.nodes[idx].left_idx  = gamma;
            self.nodes[idx].right_idx = gamma + 1;
            if(thrust::min(ij.x, ij.y) == gamma)
            {
                self.nodes[idx].left_idx += num_objects - 1;
            }
            if(thrust::max(ij.x, ij.y) == gamma + 1)
            {
                self.nodes[idx].right_idx += num_objects - 1;
            }
            self.nodes[self.nodes[idx].left_idx].parent_idx  = idx;
            self.nodes[self.nodes[idx].right_idx].parent_idx = idx;
            return;
        });
    return;
}

} // detail

template<typename Real, typename Object>
struct default_morton_code_calculator
{
    default_morton_code_calculator(aabb<Real> w): whole(w) {}
    default_morton_code_calculator()  = default;
    ~default_morton_code_calculator() = default;
    default_morton_code_calculator(default_morton_code_calculator const&) = default;
    default_morton_code_calculator(default_morton_code_calculator&&)      = default;
    default_morton_code_calculator& operator=(default_morton_code_calculator const&) = default;
    default_morton_code_calculator& operator=(default_morton_code_calculator&&)      = default;

    __device__ __host__
    inline uint32_t operator()(const Object&, const aabb<Real>& box) noexcept
    {
        auto p = centroid(box);
        p.x -= whole.lower.x;
        p.y -= whole.lower.y;
        p.z -= whole.lower.z;
        p.x /= (whole.upper.x - whole.lower.x);
        p.y /= (whole.upper.y - whole.lower.y);
        p.z /= (whole.upper.z - whole.lower.z);
        return morton_code(p);
    }
    aabb<Real> whole;
};

template<typename Real, typename Object>
using  bvh_device = detail::basic_device_bvh<Real, Object, false>;
template<typename Real, typename Object>
using cbvh_device = detail::basic_device_bvh<Real, Object, true>;

template<typename Real, typename Object, typename AABBGetter,
         typename MortonCodeCalculator = default_morton_code_calculator<Real, Object>>
class bvh
{
  public:
    using real_type   = Real;
    using index_type = std::uint32_t;
    using object_type = Object;
    using aabb_type   = aabb<real_type>;
    using node_type   = detail::node;
    using aabb_getter_type  = AABBGetter;
    using morton_code_calculator_type = MortonCodeCalculator;

  public:

    template<typename InputIterator>
    bvh(InputIterator first, InputIterator last, bool query_host_enabled = false)
        : objects_h_(first, last), objects_d_(objects_h_),
          query_host_enabled_(query_host_enabled)
    {
        this->construct();
    }

    bvh()                      = default;
    ~bvh()                     = default;
    bvh(const bvh&)            = default;
    bvh(bvh&&)                 = default;
    bvh& operator=(const bvh&) = default;
    bvh& operator=(bvh&&)      = default;

    bool  query_host_enabled() const noexcept {return query_host_enabled_;}
    bool& query_host_enabled()       noexcept {return query_host_enabled_;}

    void clear()
    {
        this->objects_h_.clear();
        this->objects_d_.clear();
        this->aabbs_h_.clear();
        this->aabbs_.clear();
        this->nodes_h_.clear();
        this->nodes_.clear();
        this->morton.clear();
        this->morton_h.clear();
        this->morton64.clear();
        this->morton64_h.clear();
        this->indices.clear();
        this->indices_h.clear();
        return ;
    }

    template<typename InputIterator>
    void assign(InputIterator first, InputIterator last)
    {
        this->objects_h_.assign(first, last);
        this->objects_d_ = this->objects_h_;
        this->construct();
        return;
    }

    bvh_device<real_type, object_type> get_device_repr()       noexcept
    {
        return bvh_device<real_type, object_type>{
            static_cast<uint32_t>(nodes_.size()),
            static_cast<uint32_t>(objects_d_.size()),
            nodes_.data().get(), aabbs_.data().get(), objects_d_.data().get()
        };
    }
    cbvh_device<real_type, object_type> get_device_repr() const noexcept
    {
        return cbvh_device<real_type, object_type>{
            static_cast<uint32_t>(nodes_.size()),
            static_cast<uint32_t>(objects_d_.size()),
            nodes_.data().get(), aabbs_.data().get(), objects_d_.data().get()
        };
    }

    void construct()
    {
        assert(objects_h_.size() == objects_d_.size());
        if(objects_h_.size() == 0u) {return;}

        const uint32_t num_objects        = objects_h_.size();
        const uint32_t num_internal_nodes = num_objects - 1;
        const uint32_t num_nodes          = num_objects * 2 - 1;

        // --------------------------------------------------------------------
        // calculate morton code of each points
        std::chrono::time_point<tempo_t> start = tempo_t::now();

        /*esta es una inicialización de el elemento identidad para hacer
        la reducción donde se calcula "whole"*/
        const auto inf = std::numeric_limits<real_type>::infinity();
        aabb_type default_aabb;
        default_aabb.upper.x = -inf; default_aabb.lower.x = inf;
        default_aabb.upper.y = -inf; default_aabb.lower.y = inf;
        default_aabb.upper.z = -inf; default_aabb.lower.z = inf;

        this->aabbs_.resize(num_nodes, default_aabb);

        /* aquí es donde se obtiene "aabb_"; se llama al operador() con cada elemento de
        objects_d_ para calcular las áreas */
        thrust::transform(this->objects_d_.begin(), this->objects_d_.end(),
                aabbs_.begin() + num_internal_nodes, aabb_getter_type());

        /* aquí se calcula "whole", que son las dimensiones de todo el área, haciendo "merge"
        de todos los puntos; esto podemos evitarlo porque en la cabecera .las tenemos max y min */
#ifdef CHECK
        aabb_type aabb_whole;

        aabb_whole.lower.x   = 715244.96;
        aabb_whole.lower.y   = 4286623.63;
        aabb_whole.lower.z   = 836.424;
        aabb_whole.upper.x   = 716057.75;
        aabb_whole.upper.y   = 4287447.70;
        aabb_whole.upper.z   = 976.790;
#else
        const aabb_type aabb_whole = thrust::reduce(
            aabbs_.begin() + num_internal_nodes, aabbs_.end(), default_aabb,
            [] __device__ (const aabb_type& lhs, const aabb_type& rhs) {
                return merge(lhs, rhs);
            });
#endif


#ifdef DEBUG
        std::cout << "aabb_whole: ";
        std::cout << aabb_whole.upper.x << "," << aabb_whole.upper.y << "," << aabb_whole.upper.z << " ";
        std::cout << aabb_whole.lower.x << "," << aabb_whole.lower.y << "," << aabb_whole.lower.z << " ";
        std::cout << std::endl;
#endif

        // thrust::device_vector<uint32_t> morton(num_objects);
        this->morton.resize(num_objects, 0u);
        /* llama al operador () de morton_code_calculator_type, el primer argumento son los
        elementos de "objects_d_" y el segundo los de "aabbs_", guardando el return en "morton" */
        thrust::transform(this->objects_d_.begin(), this->objects_d_.end(),
            aabbs_.begin() + num_internal_nodes, 
            morton.begin(),
            morton_code_calculator_type(aabb_whole));

        double mytime = cast_t(tempo_t::now() - start).count();
        std::cout << "Calculate Morton: " << mytime << " ms\n";
        double totaltime = mytime;

        // --------------------------------------------------------------------
        // sort object-indices by morton code
        start = tempo_t::now();

        // thrust::device_vector<uint32_t> indices(num_objects);
        this->indices.resize(num_objects);
        thrust::copy(thrust::make_counting_iterator<index_type>(0),
                     thrust::make_counting_iterator<index_type>(num_objects),
                     indices.begin());

#ifdef DEBUG
        for(auto item : morton){
            std::cout << item << " ";
        }
        std::cout << std::endl;

        for(auto item : indices){
            std::cout << item << " ";
        }
        std::cout << "\n";

        // for(auto itt : aabbs_){
        //     auto item = static_cast<aabb_type>(itt);
        //     std::cout << item.upper.x << "," << item.upper.y << "," << item.upper.z << " ";
        //     std::cout << item.lower.x << "," << item.lower.y << "," << item.lower.z << " ";
        //     std::cout << std::endl;
        // }
        // std::cout << std::endl;
#endif
        
        /* keep indices ascending order */
        /* ordena 3 vectores: el de los códigos Morton, el de los índices y el de BB */
        // thrust::stable_sort_by_key(
        //     morton.begin(), 
        //     morton.end(),
        //     thrust::make_zip_iterator( thrust::make_tuple(aabbs_.begin() + num_internal_nodes, indices.begin()) )
        // );

        thrust::stable_sort_by_key(thrust::device, 
            morton.begin(), 
            morton.end(),
            thrust::make_zip_iterator( thrust::make_tuple(aabbs_.begin() + num_internal_nodes, indices.begin()) ),
            // indices.begin(),
            thrust::less<uint32_t>()
        );

        
#ifdef DEBUG
        for(auto item : morton){
            std::cout << item << " ";
        }
        std::cout << std::endl;
        for(auto item : indices){
            std::cout << item << " ";
        }
        std::cout << std::endl;
#endif

        mytime = cast_t(tempo_t::now() - start).count();
        std::cout << "Sort Morton: " << mytime << " ms\n";
        totaltime += mytime;

        // --------------------------------------------------------------------
        // check morton codes are unique
        start = tempo_t::now();

        // thrust::device_vector<uint64_t> morton64(num_objects);
        this->morton64.resize(num_objects);
        /* compruebo que no hay ningun codigo morton repetido */
        const auto uniqued = thrust::unique_copy(morton.begin(), morton.end(),
                                                 morton64.begin());
        /* si no se repite ninguno, obtengo un iterador apuntando al final de "morton64" */
        const bool morton_code_is_unique = (morton64.end() == uniqued);
        /* si encuentro algún código repetido, tengo que insertar el índice al final de todos los códigos */
        if(!morton_code_is_unique)
        {
            std::cout << " Codigos repetidos -> morton 64 bits\n";
            thrust::transform(morton.begin(), morton.end(), indices.begin(),
                morton64.begin(),
                [] __device__ (const uint32_t m, const uint32_t idx)
                {
                    uint64_t m64 = m;
                    m64 <<= 32;
                    m64 |= idx;
                    return m64;
                });
        }

#ifdef DEBUG
        for(auto item : morton64){
            std::cout << item << " ";
        }
        std::cout << std::endl;
        for(auto item : indices){
            std::cout << item << " ";
        }
        std::cout << std::endl;
#endif

        mytime = cast_t(tempo_t::now() - start).count();
        std::cout << "Check Morton: " << mytime << " ms\n";
        totaltime += mytime;
        // --------------------------------------------------------------------
        // construct leaf nodes and aabbs
        start = tempo_t::now();

        node_type default_node;
        default_node.parent_idx = 0xFFFFFFFF;
        default_node.left_idx   = 0xFFFFFFFF;
        default_node.right_idx  = 0xFFFFFFFF;
        default_node.object_idx = 0xFFFFFFFF;
        this->nodes_.resize(num_nodes, default_node);

        thrust::transform(indices.begin(), indices.end(),
            this->nodes_.begin() + num_internal_nodes,
            [] __device__ (const index_type idx)
            {
                node_type n;
                n.parent_idx = 0xFFFFFFFF;
                n.left_idx   = 0xFFFFFFFF;
                n.right_idx  = 0xFFFFFFFF;
                n.object_idx = idx;
                return n;
            });
        mytime = cast_t(tempo_t::now() - start).count();
        std::cout << "Construct Leafs and AABBs: " << mytime << " ms\n";
        totaltime += mytime;

#ifdef DEBUG
        for(auto item : this->nodes_){
            auto n = static_cast<node_type>(item);
            std::cout << n.parent_idx << ", " << n.left_idx << ", ";
            std::cout << n.right_idx << ", " << n.object_idx << "\n";
        }
        std::cout << std::endl;

        for(auto itt : aabbs_){
            auto item = static_cast<aabb_type>(itt);
            std::cout << item.upper.x << "," << item.upper.y << "," << item.upper.z << " ";
            std::cout << item.lower.x << "," << item.lower.y << "," << item.lower.z << "\n";
        }
        std::cout << std::endl;
#endif


        // --------------------------------------------------------------------
        // construct internal nodes
        start = tempo_t::now();

        const auto self = this->get_device_repr();
        if(morton_code_is_unique)
        {
            const uint32_t* node_code = morton.data().get();
            detail::construct_internal_nodes(self, node_code, num_objects);
        }
        else // 64bit version
        {
            const uint64_t* node_code = morton64.data().get();
            detail::construct_internal_nodes(self, node_code, num_objects);
        }
        mytime = cast_t(tempo_t::now() - start).count();
        std::cout << "Construct Internal Nodes: " << mytime << " ms\n";
        totaltime += mytime;

#ifdef DEBUG
        for(auto item : this->nodes_){
            auto n = static_cast<node_type>(item);
            std::cout << n.parent_idx << ", " << n.left_idx << ", ";
            std::cout << n.right_idx << ", " << n.object_idx << "\n";
        }
        std::cout << std::endl;
#endif

        // --------------------------------------------------------------------
        // create AABB for each node by bottom-up strategy
        start = tempo_t::now();

        thrust::device_vector<int> flag_container(num_internal_nodes, 0);
        const auto flags = flag_container.data().get();

        thrust::for_each(thrust::device,
            thrust::make_counting_iterator<index_type>(num_internal_nodes),
            thrust::make_counting_iterator<index_type>(num_nodes),
            [self, flags] __device__ (index_type idx)
            {
                uint32_t parent = self.nodes[idx].parent_idx;
                while(parent != 0xFFFFFFFF) // means idx == 0
                {
                    const int old = atomicCAS(flags + parent, 0, 1);
                    if(old == 0)
                    {
                        // this is the first thread entered here.
                        // wait the other thread from the other child node.
                        return;
                    }
                    assert(old == 1);
                    // here, the flag has already been 1. it means that this
                    // thread is the 2nd thread. merge AABB of both childlen.

                    const auto lidx = self.nodes[parent].left_idx;
                    const auto ridx = self.nodes[parent].right_idx;
                    const auto lbox = self.aabbs[lidx];
                    const auto rbox = self.aabbs[ridx];
                    self.aabbs[parent] = merge(lbox, rbox);

                    // look the next parent...
                    parent = self.nodes[parent].parent_idx;
                }
                return;
            });
        mytime = cast_t(tempo_t::now() - start).count();
        std::cout << "Create AABB for each node: " << mytime << " ms\n";
        totaltime += mytime;
        std::cout << "CREATION takes: " << totaltime << " ms\n";

#ifdef DEBUG
        for(auto itt : aabbs_){
            auto item = static_cast<aabb_type>(itt);
            std::cout << item.upper.x << "," << item.upper.y << "," << item.upper.z << " ";
            std::cout << item.lower.x << "," << item.lower.y << "," << item.lower.z << "\n";
        }
        std::cout << std::endl;
#endif
        
        if(this->query_host_enabled_)
        {
            aabbs_h_ = aabbs_;
            nodes_h_ = nodes_;
            morton_h = morton;
            morton64_h = morton64;
            indices_h = indices;
        }
        return;
    }

    /*métodos para poder chequear los resultados en main*/
    thrust::host_vector<object_type> const& objects_host() const noexcept {return objects_h_;}
    thrust::host_vector<object_type>&       objects_host()       noexcept {return objects_h_;}
    thrust::host_vector<node_type> const& nodes_host() const noexcept {return nodes_h_;}
    thrust::host_vector<node_type>&       nodes_host()       noexcept {return nodes_h_;}
    thrust::host_vector<aabb_type> const& aabbs_host() const noexcept {return aabbs_h_;}
    thrust::host_vector<aabb_type>&       aabbs_host()       noexcept {return aabbs_h_;}
    thrust::host_vector<uint32_t> const& morton_host() const noexcept {return morton_h;}
    thrust::host_vector<uint32_t>&       morton_host()       noexcept {return morton_h;}
    thrust::host_vector<uint64_t> const& morton64_host() const noexcept {return morton64_h;}
    thrust::host_vector<uint64_t>&       morton64_host()       noexcept {return morton64_h;}
    thrust::host_vector<uint32_t> const& indices_host() const noexcept {return indices_h;}
    thrust::host_vector<uint32_t>&       indices_host()       noexcept {return indices_h;}

  private:

    thrust::host_vector  <object_type>   objects_h_;
    thrust::device_vector<object_type>   objects_d_;
    thrust::host_vector  <aabb_type>     aabbs_h_;
    thrust::device_vector<aabb_type>     aabbs_;
    thrust::host_vector  <node_type>     nodes_h_;
    thrust::device_vector<node_type>     nodes_;
    thrust::device_vector<uint32_t>     morton;
    thrust::host_vector<uint32_t>       morton_h;
    thrust::device_vector<uint64_t>     morton64;
    thrust::host_vector<uint64_t>       morton64_h;
    thrust::device_vector<uint32_t>     indices;
    thrust::host_vector<uint32_t>       indices_h;
    bool query_host_enabled_;
};

} // lbvh
#endif// LBVH_BVH_CUH
