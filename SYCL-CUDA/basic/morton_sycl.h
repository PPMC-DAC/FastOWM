#ifndef MORTON_H

#define MORTON_H

// __device__ __host__ __forceinline__
inline uint32_t SeparateBy1(uint32_t x) {
    x &= 0x0000ffffu;                  // x = ---- ---- ---- ---- fedc ba98 7654 3210
    x = (x ^ (x <<  8)) & 0x00ff00ffu; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x <<  4)) & 0x0f0f0f0fu; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x <<  2)) & 0x33333333u; // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x <<  1)) & 0x55555555u; // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    return x;
}

// __device__ __host__ __forceinline__
inline uint32_t CompactBy1(uint32_t x) {
    x &= 0x55555555;                  // x = -f-e -d-c -b-a -9-8 -7-6 -5-4 -3-2 -1-0
    x = (x ^ (x >>  1)) & 0x33333333; // x = --fe --dc --ba --98 --76 --54 --32 --10
    x = (x ^ (x >>  2)) & 0x0f0f0f0f; // x = ---- fedc ---- ba98 ---- 7654 ---- 3210
    x = (x ^ (x >>  4)) & 0x00ff00ff; // x = ---- ---- fedc ba98 ---- ---- 7654 3210
    x = (x ^ (x >>  8)) & 0x0000ffff; // x = ---- ---- ---- ---- fedc ba98 7654 3210
    return x;
}

template<typename object_t>
// __device__ __host__ __forceinline__
inline uint32_t morton2D( object_t xy , const real_t resolution = 65536.0 ) noexcept
{
  xy.x = cl::sycl::fmin(xy.x * resolution, resolution - 1.0);
  xy.y = cl::sycl::fmin(xy.y * resolution, resolution - 1.0);

  const uint32_t xx = SeparateBy1(static_cast<uint32_t>(xy.x));
  const uint32_t yy = SeparateBy1(static_cast<uint32_t>(xy.y));
  return xx * 2 + yy; // (xx << 1) | yy
}

// __device__ __host__ __forceinline__
inline uint64_t CompactBy2(uint64_t x) {
  x &= 0x1249249249249249u;                  // x = ---k --j- -i-- h--g --f- -e-- d--c --b- -a-- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
  x = (x ^ (x >> 2))  & 0x10C30C30C30C30C3u; // x = ---k ---- ji-- --hg ---- fe-- --dc ---- ba-- --98 ---- 76-- --54 ---- 32-- --10
  x = (x ^ (x >> 4))  & 0x100F00F00F00F00Fu; // x = ---k ---- ---- jihg ---- ---- fedc ---- ---- ba98 ---- ---- 7654 ---- ---- 3210
  x = (x ^ (x >> 8))  & 0x001F0000FF0000FFu; // x = ---- ---- ---k jihg ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- 7654 3210
  x = (x ^ (x >> 16)) & 0x001F00000000FFFFu; // x = ---- ---- ---k jihg ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x >> 32)) & 0x00000000001FFFFFu; // x = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---k jihg fedc ba98 7654 3210
  return x;
}

// __device__ __host__ __forceinline__
inline point_t point3D( const uint64_t code , const real_t resolution = 2097152.0 ) noexcept // 2^21
{
  const uint64_t xx = CompactBy2(code >> 2);
  const uint64_t yy = CompactBy2(code >> 1);
  const uint64_t zz = CompactBy2(code);

  point_t xyz;
  xyz.x = real_t(xx) / resolution;
  xyz.y = real_t(yy) / resolution;
  xyz.z = real_t(zz) / resolution;

  return xyz;
}

// __device__ __host__ __forceinline__
inline uint64_t SeparateBy2(uint64_t x) {
  x &= 0x00000000001FFFFFu;                  // x = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---k jihg fedc ba98 7654 3210
  x = (x ^ (x << 32)) & 0x001F00000000FFFFu; // x = ---- ---- ---k jihg ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x << 16)) & 0x001F0000FF0000FFu; // x = ---- ---- ---k jihg ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- 7654 3210
  x = (x ^ (x << 8))  & 0x100F00F00F00F00Fu; // x = ---k ---- ---- jihg ---- ---- fedc ---- ---- ba98 ---- ---- 7654 ---- ---- 3210
  x = (x ^ (x << 4))  & 0x10C30C30C30C30C3u; // x = ---k ---- ji-- --hg ---- fe-- --dc ---- ba-- --98 ---- 76-- --54 ---- 32-- --10
  x = (x ^ (x << 2))  & 0x1249249249249249u; // x = ---k --j- -i-- h--g --f- -e-- d--c --b- -a-- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
  return x;
}

template<typename object_t>
// __device__ __host__ __forceinline__
inline uint64_t morton3D( object_t xyz , const real_t resolution = 2097152.0 ) noexcept // 2^21
{
  xyz.x = cl::sycl::fmin(xyz.x * resolution, resolution - 1.0);
  xyz.y = cl::sycl::fmin(xyz.y * resolution, resolution - 1.0);
  xyz.z = cl::sycl::fmin(xyz.z * resolution, resolution - 1.0);

  const uint64_t xx = SeparateBy2(static_cast<uint64_t>(xyz.x));
  const uint64_t yy = SeparateBy2(static_cast<uint64_t>(xyz.y));
  const uint64_t zz = SeparateBy2(static_cast<uint64_t>(xyz.z));

  return xx * 4u + yy * 2u + zz; // (xx << 2) | (yy << 1) | zz
}

#endif // MORTON_H