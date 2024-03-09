#ifndef TYPES_H

#define TYPES_H

#pragma once

#include <cstdint>
#include <limits>
#include <chrono>

template<typename T>
struct vec2
{
    T x;
    T y;
};

template<typename T>
struct vec4
{
    T x;
    T y;
    T z;
    T w;
};

#ifdef USE_DOUBLE
using real_t = double;
#else
using real_t = float;
#endif

using real_v = vec2<real_t>;
using point_t = vec4<real_t>;

const auto inf = std::numeric_limits<real_t>::infinity();

struct aabb_t
{
    real_v upper;
    real_v lower;
};

struct whole_t
{
    point_t upper;
    point_t lower;
};

aabb_t default_aabb = {-inf,
                       -inf,
                       inf,
                       inf };

// using atomic_t = std::atomic<uint32_t>;

// #include <CL/sycl/atomic.hpp>
// using atomic_t = sycl::atomic<uint32_t>;

// #include <CL/sycl/atomic.hpp>
using atomic_acc = sycl::atomic_ref<uint32_t,   sycl::memory_order::relaxed,//seq_cst, 
                                                        sycl::memory_scope::device,
                                                        sycl::access::address_space::global_space>;

using atomic_wg = sycl::atomic_ref<uint32_t,   sycl::memory_order::relaxed,//seq_cst, 
                                                        sycl::memory_scope::work_group,
                                                        sycl::access::address_space::local_space>;

// using atomic_acc = sycl::atomic_ref<uint32_t,   sycl::memory_order::relaxed, 
//                                                         sycl::memory_scope::device,
//                                                         sycl::access::address_space::global_device_space>;

// using atomic_acc = sycl::atomic_ref<uint32_t,   sycl::memory_order::relaxed, 
//                                                         sycl::memory_scope::system,
//                                                         sycl::access::address_space::global_space>;

// using atomic_acc = sycl::atomic_ref<uint32_t,   sycl::memory_order::acq_rel, 
//                                                         sycl::memory_scope::device,
//                                                         sycl::access::address_space::global_device_space>;


using tempo_t = std::chrono::steady_clock;

using cast_t = std::chrono::duration<double, std::milli>;

#define TOL 1e-10f
#define FW_S32_MIN  (~0x7FFFFFFF)
#define NUM_PROCS 8

#endif //TYPES_H