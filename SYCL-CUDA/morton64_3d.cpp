#include <iostream>
#include <bitset>
#include <unistd.h>
#include <algorithm>

#define RES 2097152.0

struct double4{
	double x;
	double y;
	double z;
	double w;
};

uint64_t SeparateBy2(uint64_t x) {
  x &= 0x00000000001FFFFFu;                  // x = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---k jihg fedc ba98 7654 3210
  std::cout << std::bitset<64>(uint64_t(x)) << "\n";
  x = (x ^ (x << 32)) & 0x001F00000000FFFFu; // x = ---- ---- ---k jihg ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
  std::cout << std::bitset<64>(uint64_t(x)) << "\n";
  x = (x ^ (x << 16)) & 0x001F0000FF0000FFu; // x = ---- ---- ---k jihg ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- 7654 3210
  std::cout << std::bitset<64>(uint64_t(x)) << "\n";
  x = (x ^ (x << 8))  & 0x100F00F00F00F00Fu; // x = ---k ---- ---- jihg ---- ---- fedc ---- ---- ba98 ---- ---- 7654 ---- ---- 3210
  std::cout << std::bitset<64>(uint64_t(x)) << "\n";
  x = (x ^ (x << 4))  & 0x10C30C30C30C30C3u; // x = ---k ---- ji-- --hg ---- fe-- --dc ---- ba-- --98 ---- 76-- --54 ---- 32-- --10
  std::cout << std::bitset<64>(uint64_t(x)) << "\n";
  x = (x ^ (x << 2))  & 0x1249249249249249u; // x = ---k --j- -i-- h--g --f- -e-- d--c --b- -a-- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
  std::cout << std::bitset<64>(uint64_t(x)) << "\n\n";
  return x;
}


uint64_t CompactBy2(uint64_t x) {
  x &= 0x1249249249249249u;                  // x = ---k --j- -i-- h--g --f- -e-- d--c --b- -a-- 9--8 --7- -6-- 5--4 --3- -2-- 1--0
  x = (x ^ (x >> 2))  & 0x10C30C30C30C30C3u; // x = ---k ---- ji-- --hg ---- fe-- --dc ---- ba-- --98 ---- 76-- --54 ---- 32-- --10
  x = (x ^ (x >> 4))  & 0x100F00F00F00F00Fu; // x = ---k ---- ---- jihg ---- ---- fedc ---- ---- ba98 ---- ---- 7654 ---- ---- 3210
  x = (x ^ (x >> 8))  & 0x001F0000FF0000FFu; // x = ---- ---- ---k jihg ---- ---- ---- ---- fedc ba98 ---- ---- ---- ---- 7654 3210
  x = (x ^ (x >> 16)) & 0x001F00000000FFFFu; // x = ---- ---- ---k jihg ---- ---- ---- ---- ---- ---- ---- ---- fedc ba98 7654 3210
  x = (x ^ (x >> 32)) & 0x00000000001FFFFFu; // x = ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---k jihg fedc ba98 7654 3210
  return x;
}


uint64_t morton3D( double4 xyz , const double resolution = RES ) noexcept // 2^21
{
	xyz.x = std::min(xyz.x * resolution, resolution - 1.0);
	std::cout << std::bitset<64>(uint64_t(xyz.x)) << "\n";
	xyz.y = std::min(xyz.y * resolution, resolution - 1.0);
	std::cout << std::bitset<64>(uint64_t(xyz.y)) << "\n";
	xyz.z = std::min(xyz.z * resolution, resolution - 1.0);
	std::cout << std::bitset<64>(uint64_t(xyz.z)) << "\n\n";

  const uint64_t xx = SeparateBy2(static_cast<uint64_t>(xyz.x));
  const uint64_t yy = SeparateBy2(static_cast<uint64_t>(xyz.y));
  const uint64_t zz = SeparateBy2(static_cast<uint64_t>(xyz.z));
  return xx * 4u + yy * 2u + zz; // (xx << 2) | (yy << 1) | zz
}


int main()
{
    double4 p = {715366.194,4286696.951,844.11,0.0};
    printf("%lf %lf %lf\n", p.x,p.y,p.z);
    p.x -= 715244.96;
    p.y -= 4286623.63;
    p.z -= 836.424;
    p.x /= (716057.75 - 715244.96);
    p.y /= (4287447.70 - 4286623.63);
    p.z /= (976.790 - 836.424);
    printf("%lf %lf %lf\n", p.x,p.y,p.z);
    std::cout << uint64_t(p.x*RES) << "\n";
    std::cout << uint64_t(p.y*RES) << "\n";
    std::cout << uint64_t(p.z*RES) << "\n";

    uint64_t mortonP = morton3D(p);
    
	std::cout << std::bitset<64>(mortonP) << "\n";
    
	uint64_t dux = CompactBy2(mortonP >> 2);
    uint64_t duy = CompactBy2(mortonP >> 1);
    uint64_t duz = CompactBy2(mortonP);
    
	std::cout << std::bitset<64>(dux) << "\n";
    std::cout << std::bitset<64>(duy) << "\n";
    std::cout << std::bitset<64>(duz) << "\n";
    std::cout << dux << "\n";
    std::cout << duy << "\n";
    std::cout << duz << "\n";
    double dx = double(dux) / RES;
    double dy = double(duy) / RES;
    double dz = double(duz) / RES;
    printf("%lf %lf %lf\n", dx,dy,dz);
    dx *= (716057.75 - 715244.96);
    dy *= (4287447.70 - 4286623.63);
    dz *= (976.790 - 836.424);
    dx += 715244.96;
    dy += 4286623.63;
    dz += 836.424;
    printf("%lf %lf %lf\n", dx,dy,dz);

	return 0;
}
