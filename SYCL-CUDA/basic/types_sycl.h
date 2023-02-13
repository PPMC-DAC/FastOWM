#ifndef TYPES_H

#define TYPES_H

#pragma once

#include <cstdint>
#include <limits>
#include <chrono>

struct double2
{
    double x;
    double y;
};

struct double4
{
    double x;
    double y;
    double z;
    double w;
};

using real_t = double;
using real_v = double2;
using point_t = double4;

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
// using atomic_t = cl::sycl::atomic<uint32_t>;

// #include <CL/sycl/atomic.hpp>
using atomic_acc = sycl::atomic_ref<uint32_t,   sycl::memory_order::seq_cst, 
                                                        sycl::memory_scope::device>;
                                                        //sycl::access::address_space::global_device_space>;

using atomic_wg = sycl::atomic_ref<uint32_t,   sycl::memory_order::seq_cst, 
                                                        sycl::memory_scope::work_group>;
                                                        //, sycl::access::address_space::global_device_space>;

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

#define TOL 1e-10
#define FW_S32_MIN  (~0x7FFFFFFF)
#define NUM_PROCS 8

#endif //TYPES_H