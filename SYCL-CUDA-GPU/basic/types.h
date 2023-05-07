#ifndef TYPES_H

#define TYPES_H

#pragma once

#include <cstdint>
#include <limits>
#include <chrono>

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

using tempo_t = std::chrono::steady_clock;

using cast_t = std::chrono::duration<double, std::milli>;

#define TOL 1e-10
#define FW_S32_MIN  (~0x7FFFFFFF)
#define NUM_PROCS 8

#endif //TYPES_H